function [N1, N2, Nres] = get_OP_incNavierStokes_SUPG(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the navier stokes upwinding

% check the systems variables (u, v, p)
if sum(ismember(self.sys_variables, {'u', 'v', 'p'})) ~= 3
    error('this functions works for the field variables u, v, p')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; uv1_val = temp; up1_val = temp; uu2_val = temp; uv2_val = temp;
vu1_val = temp; vv1_val = temp; vp1_val = temp; vu2_val = temp; vv2_val = temp;
pu1_val = temp; pv1_val = temp; pp1_val = temp; pu2_val = temp; pv2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], numel(self.sys_variables));

% get the u velocities and its derivatives at the int points
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
u = self.mesh.get_fun_val(uvec, 'int'); [dudx, dudy] = self.mesh.get_1st_deriv(uvec, 'int');
[d2udx2, d2udy2] = self.mesh.get_2nd_deriv(uvec, 'int');

% get the v velocities and its derivatives at the int points
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);
v = self.mesh.get_fun_val(vvec, 'int'); [dvdx, dvdy] = self.mesh.get_1st_deriv(vvec, 'int');
[d2vdx2, d2vdy2] = self.mesh.get_2nd_deriv(vvec, 'int');

% get the pressure its derivatives at the int points
pid = find(ismember(self.sys_variables,'p')); ps = (pid-1) * n; pvec = curr_cond(:,pid);
p = self.mesh.get_fun_val(pvec, 'int'); [dpdx, dpdy] = self.mesh.get_1st_deriv(pvec, 'int');

% get the turbulence viscosity field and physical props
mue = opts.mue; rho = opts.rho; order = self.mesh.ref_elem.order;
if ~isfield(opts, 'mtvec'), opts.mtvec = zeros(1, n); end
mt = self.mesh.get_fun_val(opts.mtvec, 'int');

% loop over all the elements
count = 1;
for i=1:n_elem
    
    %% preallocate the mappings and stabilization terms
    % get the current element and its mappings
    curr_elem = self.mesh.elems(i,:);
    
    % get shortcuts to the mapping functions
    weight = self.weights(i,:)'; %inv_jac = self.inv_jac(:,:,:,i);
    dbdx = self.dbdx(:,:,i); dbdy = self.dbdy(:,:,i);
    d2bdx2 = self.d2bdx2(:,:,i); d2bdy2 = self.d2bdy2(:,:,i);
    
    % get the strong residuals in the current cell
    xres = -(mue + mt(i,:)).*(d2udx2(i,:) + d2udy2(i,:)) + rho*(u(i,:).*dudx(i,:) + v(i,:).*dudy(i,:)) + dpdx(i,:);
    yres = -(mue + mt(i,:)).*(d2vdx2(i,:) + d2vdy2(i,:)) + rho*(u(i,:).*dvdx(i,:) + v(i,:).*dvdy(i,:)) + dpdy(i,:);
    
    % calculate the support length of the element (this part may has to be reviewed)
    [tau_m, ~] = get_stabpar_wikipedia(order, rho, (mue + mt(i,:)), [u; v], weight);
    
    % loop over all of the elements nodes twice
    for k=1:n_base
        
        % define the current supg testfunction
        t_m = rho*(u(i,:).*dbdx(k,:) + v(i,:).*dbdy(k,:));
        
        for j=1:n_base
            %% calculate the u block (u test functions)
            % get the stabilization contribution (Part 1/Oseen Part)
            s_uu1 = tau_m.*(-(mue + mt(i,:)).*(d2bdx2(j,:) + d2bdy2(j,:)) + rho*(u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
            s_up1 = tau_m.*(dbdx(j,:)).*t_m; % pressure part of the momentum residual
            
            % get the stabilization contribution (Part 2)
            s_uu2 = tau_m.*(rho*(b(j,:).*dudx(i,:)).*t_m + rho*(xres).*b(j,:).*dbdx(k,:));
            s_uv2 = tau_m.*(rho*(b(j,:).*dudy(i,:)).*t_m + rho*(xres).*b(j,:).*dbdy(k,:));
            
            % integrate the residuals function
            uu1_val(count) = s_uu1*weight; uu2_val(count) = s_uu2*weight;
            uv2_val(count) = s_uv2*weight;
            up1_val(count) = s_up1*weight;
            
            %% calculate the v block (v test functions)
            % get the stabilization contribution (Part 1/ Oseen Part)
            s_vv1 = tau_m.*(-(mue + mt(i,:)).*(d2bdx2(j,:) + d2bdy2(j,:)) + rho*(u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
            s_vp1 = tau_m.*(dbdy(j,:)).*t_m; % pressure part of the momentum residual
            
            % get the stabilization contribution (Part 2)
            s_vu2 = tau_m.*(rho*(b(j,:).*dvdx(i,:)).*t_m + rho*(yres).*b(j,:).*dbdx(k,:));
            s_vv2 = tau_m.*(rho*(b(j,:).*dvdy(i,:)).*t_m + rho*(yres).*b(j,:).*dbdy(k,:));
            
            % integrate the residuals function
            vu2_val(count) = s_vu2*weight;
            vv1_val(count) = s_vv1*weight; vv2_val(count) = s_vv2*weight;
            vp1_val(count) = s_vp1*weight;
            
            %% calculate the p block (p test functions)
            % get the stabilization contribution (Part 1/ Oseen Part)
            s_pu1 = tau_m.*(-(mue + mt(i,:)).*(d2bdx2(j,:) + d2bdy2(j,:)) + rho*(u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:))).*dbdx(k,:); % 
            s_pv1 = tau_m.*(-(mue + mt(i,:)).*(d2bdx2(j,:) + d2bdy2(j,:)) + rho*(u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:))).*dbdy(k,:); % 
            s_pp1 = tau_m.*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:));
            
            % get the stabilization contribution (Part 2)
            s_pu2 = tau_m.*(rho*(b(j,:).*dudx(i,:).*dbdx(k,:) + b(j,:).*dvdx(i,:).*dbdy(k,:)));
            s_pv2 = tau_m.*(rho*(b(j,:).*dudy(i,:).*dbdx(k,:) + b(j,:).*dvdy(i,:).*dbdy(k,:)));
            
            % integrate the residuals function
            pu1_val(count) = s_pu1*weight; pu2_val(count) = s_pu2*weight;
            pv1_val(count) = s_pv1*weight; pv2_val(count) = s_pv2*weight;
            pp1_val(count) = s_pp1*weight;
            
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