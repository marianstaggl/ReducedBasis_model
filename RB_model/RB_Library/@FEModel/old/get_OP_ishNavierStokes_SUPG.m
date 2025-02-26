function [N, Nres] = get_OP_ishNavierStokes_SUPG(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the navier stokes upwinding

% check the systems variables (u, v, p)
if sum(ismember(self.sys_variables, {'u', 'v', 'h'})) ~= 3
    error('this functions works for the field variables u, v, and h')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; order = self.mesh.ref_elem.order; gamma = 1.4;
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; uv1_val = temp; uh1_val = temp; uu2_val = temp; uv2_val = temp; uh2_val = temp;
vu1_val = temp; vv1_val = temp; vh1_val = temp; vu2_val = temp; vv2_val = temp; vh2_val = temp;
hu1_val = temp; hv1_val = temp; hh1_val = temp; hu2_val = temp; hv2_val = temp; hh2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], 3);

% get the starting points of each block and sort them into vectors
hid = find(ismember(self.sys_variables,'h')); hs = (hid-1) * n; hvec = curr_cond(:,hid);
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);

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
    uval = uvec(curr_elem); vval = vvec(curr_elem); hval = hvec(curr_elem);
    
    % get the current velocitys at the quad points
    u = sum(b.*uval,1); dudx = sum(dbdx.*uval,1); dudy = sum(dbdy.*uval,1); 
    d2udx2 = sum(d2bdx2.*uval,1); d2udy2 = sum(d2bdy2.*uval,1);
    v = sum(b.*vval,1); dvdx = sum(dbdx.*vval,1); dvdy = sum(dbdy.*vval,1); 
    d2vdx2 = sum(d2bdx2.*vval,1); d2vdy2 = sum(d2bdy2.*vval,1);
    h = sum(b.*hval,1); dhdx = sum(dbdx.*hval,1); dhdy = sum(dbdy.*hval,1);
    
    % get the strong residuals in the current cell
    xres = -mue*(d2udx2 + d2udy2) + u.*dudx + v.*dudy + dhdx;
    yres = -mue*(d2vdx2 + d2vdy2) + u.*dvdx + v.*dvdy + dhdy;
    cres = u.*dhdx + v.*dhdy + (gamma-1)*h.*(dudx + dvdy);
    
    % calculate the support length of the element (this part may has to be reviewed)
    [tau_m, tau_c] = get_stabpar_wikipedia(order, 1, mue, [u; v], weight);
    %[tau_m, tau_c] = get_stabpar_redbkit(order, 1e-5, mue, [u; v], inv_jac);
    % tau_m = 1e-2; tau_c = 1e-4; %1e-8; %0.5; tau_c = 2; %1/diam;
    
    % loop over all of the elements nodes twice
    for k=1:n_base
        
        % define the current supg testfunction
        t_m = u.*dbdx(k,:) + v.*dbdy(k,:);
        c_mx = ((gamma-1)*h.*dbdx(k,:)); c_my = ((gamma-1)*h.*dbdy(k,:)); c_mc = (u.*dbdx(k,:) + v.*dbdy(k,:));
        
        for j=1:n_base
            
            %% calculate the u block (u test functions)
            % get the grad div stabilization (Part 1/continuity residual)
            g_uu1 = tau_c.*((gamma-1)*h.*dbdx(j,:)).*c_mx; %dbdx(j,:).*dbdx(k,:));
            g_uv1 = tau_c.*((gamma-1)*h.*dbdy(j,:)).*c_mx; %dbdy(j,:).*dbdx(k,:));
            g_uh1 = tau_c.*(u.*dbdx(j,:) + v.*dbdy(j,:)).*c_mx;
            
            % get the grad div stabilization (Part 2/continuity residual)
            g_uu2 = tau_c.*(b(j,:).*dhdx).*c_mx;
            g_uv2 = tau_c.*(b(j,:).*dhdy).*c_mx;
            g_hh2 = tau_c.*(((gamma-1)*b(j,:).*(dudx + dvdy)).*c_mx + ((gamma-1)*b(j,:).*(cres).*dbdx(k,:)));
            
            % get the stabilization contribution (Part 1/Oseen Part)
            s_uu1 = tau_m.*(-mue*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
            s_uh1 = tau_m.*(dbdx(j,:)).*t_m; % pressure part of the momentum residual
            
            % get the stabilization contribution (Part 2)
            s_uu2 = tau_m.*((b(j,:).*dudx).*t_m + (xres).*b(j,:).*dbdx(k,:));
            s_uv2 = tau_m.*((b(j,:).*dudy).*t_m + (xres).*b(j,:).*dbdy(k,:));
            
            % integrate the residuals function
            uu1_val(count) = sum(weight.*(s_uu1 + g_uu1)); uu2_val(count) = sum(weight.*(s_uu2 + g_uu2));
            uv1_val(count) = sum(weight.*(g_uv1));         uv2_val(count) = sum(weight.*(s_uv2 + g_uv2));
            uh1_val(count) = sum(weight.*(s_uh1 + g_uh1)); uh2_val(count) = sum(weight.*(g_hh2));
            
            %% calculate the v block (v test functions)
            % get the grad div stabilization (continuity residual)
            g_vu1 = tau_c.*((gamma-1)*h.*dbdx(j,:)).*c_my; %dbdx(j,:).*dbdy(k,:));
            g_vv1 = tau_c.*((gamma-1)*h.*dbdy(j,:)).*c_my; %dbdy(j,:).*dbdy(k,:));
            g_vh1 = tau_c.*(u.*dbdx(j,:) + v.*dbdy(j,:)).*c_my;
            
            % get the grad div stabilization (Part 2/continuity residual)
            g_vu2 = tau_c.*(b(j,:).*dhdx).*c_my;
            g_vv2 = tau_c.*(b(j,:).*dhdy).*c_my;
            g_vh2 = tau_c.*(((gamma-1)*b(j,:).*(dudx + dvdy)).*c_my + ((gamma-1)*b(j,:).*(cres).*dbdy(k,:)));
            
            % get the stabilization contribution (Part 1/ Oseen Part)
            s_vv1 = tau_m.*(-mue*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
            s_vh1 = tau_m.*(dbdy(j,:)).*t_m; % pressure part of the momentum residual
            
            % get the stabilization contribution (Part 2)
            s_vu2 = tau_m.*((b(j,:).*dvdx).*t_m + (yres).*b(j,:).*dbdx(k,:));
            s_vv2 = tau_m.*((b(j,:).*dvdy).*t_m + (yres).*b(j,:).*dbdy(k,:));
            
            % integrate the residuals function
            vu1_val(count) = sum(weight.*(g_vu1));         vu2_val(count) = sum(weight.*(s_vu2 + g_vu2)); 
            vv1_val(count) = sum(weight.*(s_vv1 + g_vv1)); vv2_val(count) = sum(weight.*(s_vv2 + g_vv2));
            vh1_val(count) = sum(weight.*(s_vh1 + g_vh1)); vh2_val(count) = sum(weight.*(g_vh2));
            
            %% calculate the p block (p test functions)
            % get the stabilization contribution (Part 1/continuity residual
            g_hu1 = tau_c.*((gamma-1)*h.*dbdx(j,:)).*c_mc;
            g_hv1 = tau_c.*((gamma-1)*h.*dbdy(j,:)).*c_mc;
            g_hh1 = tau_c.*(u.*dbdx(j,:) + v.*dbdy(j,:)).*c_mc;
            
            % get the grad div stabilization (Part 2/continuity residual)
            g_hu2 = tau_c.*((b(j,:).*dhdx).*c_mc + (b(j,:).*(cres)).*dbdx(k,:));
            g_hv2 = tau_c.*((b(j,:).*dhdy).*c_mc + (b(j,:).*(cres)).*dbdy(k,:));
            g_hh2 = tau_c.*((gamma-1)*b(j,:).*(dudx + dvdy)).*c_mc;
            
            % get the stabilization contribution (Part 1/ Oseen Part)
            s_hu1 = tau_m.*(-mue*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*dbdx(k,:); % 
            s_hv1 = tau_m.*(-mue*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*dbdy(k,:); % 
            s_hh1 = tau_m.*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:));
            
            % get the stabilization contribution (Part 2)
            s_hu2 = tau_m.*((b(j,:).*dudx.*dbdx(k,:) + b(j,:).*dvdx.*dbdy(k,:)));
            s_hv2 = tau_m.*((b(j,:).*dudy.*dbdx(k,:) + b(j,:).*dvdy.*dbdy(k,:)));
            
            % integrate the residuals function
            hu1_val(count) = sum(weight.*(s_hu1 + g_hu1)); hu2_val(count) = sum(weight.*(s_hu2 + g_hu2));
            hv1_val(count) = sum(weight.*(s_hv1 + g_hv1)); hv2_val(count) = sum(weight.*(s_hv2 + g_hv2));
            hh1_val(count) = sum(weight.*(s_hh1 + g_hh1)); hh2_val(count) = sum(weight.*(g_hh2));
            
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
v1_glob = [uu1_val,    uv1_val,    uh1_val,    vu1_val,    vv1_val,    vh1_val,    hu1_val,    hv1_val,    hh1_val   ];
k1_glob = [k_idx + us, k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + vs, k_idx + hs, k_idx + hs, k_idx + hs];
j1_glob = [j_idx + us, j_idx + vs, j_idx + hs, j_idx + us, j_idx + vs, j_idx + hs, j_idx + us, j_idx + vs, j_idx + hs];

% reshape the vectors for the sparse matrix assembly
v2_glob = [uu2_val,    uv2_val,    vu2_val,    vv2_val,    hu2_val,    hv2_val   ];
k2_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + hs, k_idx + hs];
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
t_m = ((rho*vecnorm(vel)/h).^2 + ck*(mue/(h^2))^2).^(-1/2);

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