function [N, Nres] = get_OP_incNavierStokes_turb_SUPG(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the navier stokes upwinding

% check the systems variables (u, v, p)
if sum(ismember(self.sys_variables, {'u', 'v', 'p', 'nt'})) ~= 4
    error('this functions works for the field variables u, v, p and nt')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; rho = opts.rho; order = self.mesh.ref_elem.order;
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the residuals vector
% N2res = zeros(n*numel(self.sys_variables),1);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; uv1_val = temp; up1_val = temp; uu2_val = temp; uv2_val = temp;
vu1_val = temp; vv1_val = temp; vp1_val = temp; vu2_val = temp; vv2_val = temp;
pu1_val = temp; pv1_val = temp; pp1_val = temp; pu2_val = temp; pv2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], 4);

% get the starting points of each block and sort them into vectors
pid = find(ismember(self.sys_variables,'p')); ps = (pid-1) * n; pvec = curr_cond(:,pid);
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);

% calculate the mue_t vector from nue_t and density
ntid = ismember(self.sys_variables, 'nt'); ntvec = curr_cond(:,ntid);
cv1 = 7.1; X  = ntvec./(mue./rho); fv1 = X.^3./(X.^3 - cv1^3); 
mtvec = rho.*ntvec.*fv1;

% loop over all the elements
count = 1;
for i=1:n_elem
    
    %% preallocate the mappings and stabilization terms
    % get the current element and its mappings
    curr_elem = self.mesh.elems(i,:);
    
    % get shortcuts to the mapping functions
    weight = self.weights(i,:); %inv_jac = self.inv_jac(:,:,:,i);
    dbdx = self.dbdx(:,:,i); dbdy = self.dbdy(:,:,i);
    d2bdx2 = self.d2bdx2(:,:,i); d2bdy2 = self.d2bdy2(:,:,i);
    uval = uvec(curr_elem); vval = vvec(curr_elem); 
    pval = pvec(curr_elem); mtval = mtvec(curr_elem);
    
    % get the current velocitys at the quad points
    u = sum(b.*uval,1); dudx = sum(dbdx.*uval,1); dudy = sum(dbdy.*uval,1); 
    d2udx2 = sum(d2bdx2.*uval,1); d2udy2 = sum(d2bdy2.*uval,1);
    v = sum(b.*vval,1); dvdx = sum(dbdx.*vval,1); dvdy = sum(dbdy.*vval,1); 
    d2vdx2 = sum(d2bdx2.*vval,1); d2vdy2 = sum(d2bdy2.*vval,1);
    p = sum(b.*pval,1); dpdx = sum(dbdx.*pval,1); dpdy = sum(dbdy.*pval,1);
    mt = sum(b.*mtval,1);
    
    % get the strong residuals in the current cell
    xres = -(mue + mt).*(d2udx2 + d2udy2) + rho*(u.*dudx + v.*dudy) + dpdx;
    yres = -(mue + mt).*(d2vdx2 + d2vdy2) + rho*(u.*dvdx + v.*dvdy) + dpdy;
    
    % calculate the support length of the element (this part may has to be reviewed)
    [tau_m, ~] = get_stabpar_wikipedia(order, rho, mue + mt, [u; v], weight);
%     [tau_m, tau_c] = get_stabpar_redbkit(order, rho, mue, [u; v], inv_jac);
    
    % loop over all of the elements nodes twice
    for k=1:n_base
        
        % define the current supg testfunction
        t_m = rho*(u.*dbdx(k,:) + v.*dbdy(k,:));
        
        for j=1:n_base
            %% calculate the u block (u test functions)
            % get the stabilization contribution (Part 1/Oseen Part)
            s_uu1 = tau_m.*(-(mue + mt).*(d2bdx2(j,:) + d2bdy2(j,:)) + rho*(u.*dbdx(j,:) + v.*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
            s_up1 = tau_m.*(dbdx(j,:)).*t_m; % pressure part of the momentum residual
            
            % get the stabilization contribution (Part 2)
            s_uu2 = tau_m.*(rho*(b(j,:).*dudx).*t_m + rho*(xres).*b(j,:).*dbdx(k,:));
            s_uv2 = tau_m.*(rho*(b(j,:).*dudy).*t_m + rho*(xres).*b(j,:).*dbdy(k,:));
            
            % integrate the residuals function
            uu1_val(count) = sum(weight.*(s_uu1)); uu2_val(count) = sum(weight.*(s_uu2));
            uv2_val(count) = sum(weight.*(s_uv2));
            up1_val(count) = sum(weight.*(s_up1));
            
            %% calculate the v block (v test functions)
            % get the stabilization contribution (Part 1/ Oseen Part)
            s_vv1 = tau_m.*(-(mue + mt).*(d2bdx2(j,:) + d2bdy2(j,:)) + rho*(u.*dbdx(j,:) + v.*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
            s_vp1 = tau_m.*(dbdy(j,:)).*t_m; % pressure part of the momentum residual
            
            % get the stabilization contribution (Part 2)
            s_vu2 = tau_m.*(rho*(b(j,:).*dvdx).*t_m + rho*(yres).*b(j,:).*dbdx(k,:));
            s_vv2 = tau_m.*(rho*(b(j,:).*dvdy).*t_m + rho*(yres).*b(j,:).*dbdy(k,:));
            
            % integrate the residuals function
            vu2_val(count) = sum(weight.*(s_vu2)); 
            vv1_val(count) = sum(weight.*(s_vv1)); vv2_val(count) = sum(weight.*(s_vv2));
            vp1_val(count) = sum(weight.*(s_vp1));
            
            %% calculate the p block (p test functions)
            % get the stabilization contribution (Part 1/ Oseen Part)
            s_pu1 = tau_m.*(-(mue + mt).*(d2bdx2(j,:) + d2bdy2(j,:)) + rho*(u.*dbdx(j,:) + v.*dbdy(j,:))).*dbdx(k,:); % 
            s_pv1 = tau_m.*(-(mue + mt).*(d2bdx2(j,:) + d2bdy2(j,:)) + rho*(u.*dbdx(j,:) + v.*dbdy(j,:))).*dbdy(k,:); % 
            s_pp1 = tau_m.*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:));
            
            % get the stabilization contribution (Part 2)
            s_pu2 = tau_m.*(rho*(b(j,:).*dudx.*dbdx(k,:) + b(j,:).*dvdx.*dbdy(k,:)));
            s_pv2 = tau_m.*(rho*(b(j,:).*dudy.*dbdx(k,:) + b(j,:).*dvdy.*dbdy(k,:)));
            
            % integrate the residuals function
            pu1_val(count) = sum(weight.*(s_pu1)); pu2_val(count) = sum(weight.*(s_pu2));
            pv1_val(count) = sum(weight.*(s_pv1)); pv2_val(count) = sum(weight.*(s_pv2));
            pp1_val(count) = sum(weight.*(s_pp1));
            
            %% do other stuff
            % get the global indices for the assembly
            j_idx(count) = curr_elem(j); k_idx(count) = curr_elem(k);
            
            % count one up
            count = count + 1;
        end
    end
end

%% form a sparse matrix from the galerkin terms
% reshape the vectors for the sparse matrix assembly
v1_glob = [uu1_val,    uv1_val,    up1_val,    vu1_val,    vv1_val,    vp1_val,    pu1_val,    pv1_val,    pp1_val   ];
k1_glob = [k_idx + us, k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + vs, k_idx + ps, k_idx + ps, k_idx + ps];
j1_glob = [j_idx + us, j_idx + vs, j_idx + ps, j_idx + us, j_idx + vs, j_idx + ps, j_idx + us, j_idx + vs, j_idx + ps];

% reshape the vectors for the sparse matrix assembly
v2_glob = [uu2_val,    uv2_val,    vu2_val,    vv2_val,    pu2_val,    pv2_val   ];
k2_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + ps, k_idx + ps];
j2_glob = [j_idx + us, j_idx + vs, j_idx + us, j_idx + vs, j_idx + us, j_idx + vs];

% get a sparse global matrix
N1 = sparse(k1_glob, j1_glob, v1_glob, numel(curr_cond), numel(curr_cond));
N2 = sparse(k2_glob, j2_glob, v2_glob, numel(curr_cond), numel(curr_cond));

% get the residuals at the nodes
Nres = N1*reshape(curr_cond,[],1);
N = N1 + opts.lin*N2;

end

function [t_m, t_c] = get_stabpar_redbkit(order, rho, mue, vel, inv_jac)
% get the stabilization parameters for the current cell (Version: Redbkit
% documentary). G and g are metric Tensors of the current cell and
% represent a measure for the cell "size"

% preallocate the arrays
t_m = zeros(1,size(inv_jac,3)); t_c = zeros(1,size(inv_jac,3));

% loop over every integration point
for i=1:size(inv_jac,3)
    
    % get the inverse jacobian of the current int point
    c_invj = inv_jac(:,:,i);
    
    % get the metric tensor g and G (source: redkit github docu)
    g = sum(c_invj,2); G = c_invj * c_invj';
    
    % define the constant C
    c = 60*2^(order-2);
    
    % calculate the parameter t_m
    t_m(i) = ((rho^2)*dot(vel(:,i), G*vel(:,i)) + c*(mue^2)*sum(G.*G,'all')).^(-1/2);
    
    % calculate the parameter t_c
    t_c(i) = (t_m(i)*dot(g,g)).^(-1);
end
end

function [t_m, t_c] = get_stabpar_wikipedia(order, rho, mue, vel, weights)
% get the stabilization parameters for the current cell (Version:
% Wikipedia).

% calculate the "length" of the element (diameter of a circle)
h = 2*sqrt(sum(weights)/pi);

% get the constant ck
ck = 60*2^(order-2);

% get the parameter tau_m
t_m = ((rho*vecnorm(vel)/h).^2 + ck*(mue/(h^2)).^2).^(-1/2);

% calculate the parameter t_c
t_c = (h^2)./t_m;

end


% 
% % get the weighted residuals
% c_xres = tau_m.*(-mue*(d2udx2 + d2udy2) + rho*(u.*dudx + v.*dudy) + dpdx ).*t_m + tau_c.*( dudx + dvdy ).*dbdx(k,:);
% c_yres = tau_m.*(-mue*(d2vdx2 + d2vdy2) + rho*(u.*dvdx + v.*dvdy) + dpdy ).*t_m + tau_c.*( dudx + dvdy ).*dbdy(k,:);
% c_cres = tau_m.*((-mue*(d2udx2 + d2udy2) + rho*(u.*dudx + v.*dudy) + dpdx ).*dbdx(k,:) + (-mue*(d2vdx2 + d2vdy2) + rho*(u.*dvdx + v.*dvdy) + dpdy ).*dbdy(k,:));
% 
% % integrate the current residual
% N2res(curr_elem(k) + us) = N2res(curr_elem(k) + us) + sum(c_xres.*weight);
% N2res(curr_elem(k) + vs) = N2res(curr_elem(k) + vs) + sum(c_yres.*weight);
% N2res(curr_elem(k) + ps) = N2res(curr_elem(k) + ps) + sum(c_cres.*weight);

%     % get the current velocitys and their derivatives at the quad points
%     u = (b'*uvec(curr_elem))'; v = (b'*vvec(curr_elem))';
%     dpdx = (dbdx'*pvec(curr_elem))'; dpdy = (dbdy'*pvec(curr_elem))';
%     dudx = (dbdx'*uvec(curr_elem))'; dudy = (dbdy'*uvec(curr_elem))';
%     dvdx = (dbdx'*vvec(curr_elem))'; dvdy = (dbdy'*vvec(curr_elem))';
%     d2udx2 = (d2bdx2'*uvec(curr_elem))'; d2udy2 = (d2bdy2'*uvec(curr_elem))';
%     d2vdx2 = (d2bdx2'*vvec(curr_elem))'; d2vdy2 = (d2bdy2'*vvec(curr_elem))';