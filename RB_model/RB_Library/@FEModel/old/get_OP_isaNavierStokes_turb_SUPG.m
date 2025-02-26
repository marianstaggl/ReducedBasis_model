function [N, Nres] = get_OP_isaNavierStokes_turb_SUPG(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the navier stokes upwinding

% check the systems variables (u, v, a, nt)
if sum(ismember(self.sys_variables, {'u', 'v', 'a', 'nt'})) ~= 4
    error('this functions works for the field variables u, v, a and nt')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; order = self.mesh.ref_elem.order; gamma = 1.4;
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; ua1_val = temp; uu2_val = temp; uv2_val = temp; ua2_val = temp;
vv1_val = temp; va1_val = temp; vu2_val = temp; vv2_val = temp; va2_val = temp;
au1_val = temp; av1_val = temp; aa1_val = temp; au2_val = temp; av2_val = temp; aa2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], 4);

% get the starting points of each block and sort them into vectors
aid = find(ismember(self.sys_variables,'a')); as = (aid-1) * n; avec = curr_cond(:,aid);
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);

% calculate the mue_t vector from nue_t and density
ntid = ismember(self.sys_variables, 'nt'); ntvec = curr_cond(:,ntid); 
rho_ref = 1.1841; p_ref = 101325; cv1 = 7.1;
rho = (avec.^2/(gamma*p_ref)*(rho_ref.^(gamma))).^(1/(gamma-1));
X  = ntvec./(mue./rho); fv1 = X.^3./(X.^3 - cv1^3); mtvec = rho.*ntvec.*fv1;

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
    aval = avec(curr_elem); mtval = mtvec(curr_elem);
    
    % get the current velocitys at the quad points
    u = sum(b.*uval,1); dudx = sum(dbdx.*uval,1); dudy = sum(dbdy.*uval,1); 
    d2udx2 = sum(d2bdx2.*uval,1); d2udy2 = sum(d2bdy2.*uval,1);
    v = sum(b.*vval,1); dvdx = sum(dbdx.*vval,1); dvdy = sum(dbdy.*vval,1); 
    d2vdx2 = sum(d2bdx2.*vval,1); d2vdy2 = sum(d2bdy2.*vval,1);
    a = sum(b.*aval,1); dadx = sum(dbdx.*aval,1); dady = sum(dbdy.*aval,1);
    mt = sum(b.*mtval,1);
    
    % get the strong residuals in the current cell
    xres = -(mue + mt).*(d2udx2 + d2udy2) + u.*dudx + v.*dudy + (2/(gamma-1))*a.*dadx;
    yres = -(mue + mt).*(d2vdx2 + d2vdy2) + u.*dvdx + v.*dvdy + (2/(gamma-1))*a.*dady;
    
    % calculate the support length of the element (this part may has to be reviewed)
    [tau_m, ~] = get_stabpar_wikipedia(order, 1, mue, [u; v], weight);
    %[tau_m, tau_c] = get_stabpar_redbkit(order, 1e-5, mue, [u; v], inv_jac);
    
    % loop over all of the elements nodes twice
    for k=1:n_base
        
        % define the current supg testfunction
        t_m = u.*dbdx(k,:) + v.*dbdy(k,:);
        c_mx = ((2/(gamma-1))*a.*dbdx(k,:));
        c_my = ((2/(gamma-1))*a.*dbdy(k,:));
        
        for j=1:n_base
            
            %% calculate the u block (u test functions)
            % get the stabilization contribution (Part 1/Oseen Part)
            s_uu1 = tau_m.*(-(mue + mt).*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
            s_uh1 = tau_m.*(2/(gamma-1)*a.*dbdx(j,:)).*t_m; % pressure part of the momentum residual
            
            % get the stabilization contribution (Part 2)
            s_uu2 = tau_m.*((b(j,:).*dudx).*t_m + (xres).*b(j,:).*dbdx(k,:));
            s_uv2 = tau_m.*((b(j,:).*dudy).*t_m + (xres).*b(j,:).*dbdy(k,:));
            s_ua2 = tau_m.*(2/(gamma-1)*b(j,:).*dadx).*t_m;
            
            % integrate the residuals function
            uu1_val(count) = sum(weight.*(s_uu1)); uu2_val(count) = sum(weight.*(s_uu2));
            uv2_val(count) = sum(weight.*(s_uv2));
            ua1_val(count) = sum(weight.*(s_uh1)); ua2_val(count) = sum(weight.*(s_ua2));
            
            %% calculate the v block (v test functions)
            % get the stabilization contribution (Part 1/ Oseen Part)
            s_vv1 = tau_m.*(-(mue + mt).*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
            s_va1 = tau_m.*(2/(gamma-1)*a.*dbdy(j,:)).*t_m; % pressure part of the momentum residual
            
            % get the stabilization contribution (Part 2)
            s_vu2 = tau_m.*((b(j,:).*dvdx).*t_m + (yres).*b(j,:).*dbdx(k,:));
            s_vv2 = tau_m.*((b(j,:).*dvdy).*t_m + (yres).*b(j,:).*dbdy(k,:));
            s_va2 = tau_m.*(2/(gamma-1)*b(j,:).*dady).*t_m;
            
            % integrate the residuals function
            vu2_val(count) = sum(weight.*(s_vu2)); 
            vv1_val(count) = sum(weight.*(s_vv1)); vv2_val(count) = sum(weight.*(s_vv2));
            va1_val(count) = sum(weight.*(s_va1)); va2_val(count) = sum(weight.*(s_va2));
            
            %% calculate the p block (p test functions)
            % get the stabilization contribution (Part 1/ Oseen Part)
            s_au1 = tau_m.*(-(mue + mt).*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*c_mx; % dbdx(k,:)
            s_av1 = tau_m.*(-(mue + mt).*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*c_my; % dbdy(k,:)
            s_aa1 = tau_m.*(2/(gamma-1)*a.*(dbdx(j,:).*c_mx + dbdy(j,:).*c_my)); % tau_m.*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:));
            
            % get the stabilization contribution (Part 2)
            s_au2 = tau_m.*(b(j,:).*dudx.*c_mx);
            s_av2 = tau_m.*(b(j,:).*dvdy.*c_my);
            s_aa2 = tau_m.*(2/(gamma-1)*b(j,:).*(dadx.*c_mx + dady.*c_my));
            
            % integrate the residuals function
            au1_val(count) = sum(weight.*(s_au1)); au2_val(count) = sum(weight.*(s_au2));
            av1_val(count) = sum(weight.*(s_av1)); av2_val(count) = sum(weight.*(s_av2));
            aa1_val(count) = sum(weight.*(s_aa1)); aa2_val(count) = sum(weight.*(s_aa2));
            
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
v1_glob = [uu1_val,    ua1_val,    vv1_val,    va1_val,    au1_val,    av1_val,    aa1_val   ];
k1_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + as, k_idx + as, k_idx + as];
j1_glob = [j_idx + us, j_idx + as, j_idx + vs, j_idx + as, j_idx + us, j_idx + vs, j_idx + as];

% reshape the vectors for the sparse matrix assembly
v2_glob = [uu2_val,    uv2_val,    vu2_val,    vv2_val,    au2_val,    av2_val   , aa2_val   ];
k2_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + as, k_idx + as, k_idx + as];
j2_glob = [j_idx + us, j_idx + vs, j_idx + us, j_idx + vs, j_idx + us, j_idx + vs, j_idx + as];

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

%{
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

        c_mx = ((gamma-1)/2*a.*dbdx(k,:)); 
        c_my = ((gamma-1)/2*a.*dbdy(k,:)); 
        c_mc = (u.*dbdx(k,:) + v.*dbdy(k,:));


%             calculate the u block (u test functions)
%             % get the grad div stabilization (Part 1/continuity residual)
%             g_uu1 = tau_c.*((gamma-1)/2*a.*dbdx(j,:)).*c_mx; %dbdx(j,:).*dbdx(k,:));
%             g_uv1 = tau_c.*((gamma-1)/2*a.*dbdy(j,:)).*c_mx; %dbdy(j,:).*dbdx(k,:));
%             g_ua1 = tau_c.*(u.*dbdx(j,:) + v.*dbdy(j,:)).*c_mx;
%             
%             % get the grad div stabilization (Part 2/continuity residual)
%             g_uu2 = tau_c.*(b(j,:).*dadx).*c_mx;
%             g_uv2 = tau_c.*(b(j,:).*dady).*c_mx;
%             g_aa2 = tau_c.*(((gamma-1)/2*b(j,:).*(dudx + dvdy)).*c_mx + ((gamma-1)/2*b(j,:).*(cres).*dbdx(k,:)));
%             
%             % get the stabilization contribution (Part 1/Oseen Part)
%             s_uu1 = tau_m.*(-mue*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
%             s_uh1 = tau_m.*((gamma-1)/2*a.*dbdx(j,:)).*t_m; % pressure part of the momentum residual
%             
%             % get the stabilization contribution (Part 2)
%             s_uu2 = tau_m.*((b(j,:).*dudx).*t_m + (xres).*b(j,:).*dbdx(k,:));
%             s_uv2 = tau_m.*((b(j,:).*dudy).*t_m + (xres).*b(j,:).*dbdy(k,:));
%             s_ua2 = tau_m.*((gamma-1)/2*b(j,:).*dadx).*t_m;
%             
%             % integrate the residuals function
%             uu1_val(count) = sum(weight.*(s_uu1 + g_uu1)); uu2_val(count) = sum(weight.*(s_uu2 + g_uu2));
%             uv1_val(count) = sum(weight.*(g_uv1));         uv2_val(count) = sum(weight.*(s_uv2 + g_uv2));
%             ua1_val(count) = sum(weight.*(s_uh1 + g_ua1)); ua2_val(count) = sum(weight.*(s_ua2 + g_aa2));
%             
%             %% calculate the v block (v test functions)
%             % get the grad div stabilization (continuity residual)
%             g_vu1 = tau_c.*((gamma-1)/2*a.*dbdx(j,:)).*c_my; %dbdx(j,:).*dbdy(k,:));
%             g_vv1 = tau_c.*((gamma-1)/2*a.*dbdy(j,:)).*c_my; %dbdy(j,:).*dbdy(k,:));
%             g_va1 = tau_c.*(u.*dbdx(j,:) + v.*dbdy(j,:)).*c_my;
%             
%             % get the grad div stabilization (Part 2/continuity residual)
%             g_vu2 = tau_c.*(b(j,:).*dadx).*c_my;
%             g_vv2 = tau_c.*(b(j,:).*dady).*c_my;
%             g_va2 = tau_c.*(((gamma-1)/2*b(j,:).*(dudx + dvdy)).*c_my + ((gamma-1)/2*b(j,:).*(cres).*dbdy(k,:)));
%             
%             % get the stabilization contribution (Part 1/ Oseen Part)
%             s_vv1 = tau_m.*(-mue*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*t_m; % convective and diffusive part of the momentum residual 
%             s_va1 = tau_m.*((gamma-1)/2*a.*dbdy(j,:)).*t_m; % pressure part of the momentum residual
%             
%             % get the stabilization contribution (Part 2)
%             s_vu2 = tau_m.*((b(j,:).*dvdx).*t_m + (yres).*b(j,:).*dbdx(k,:));
%             s_vv2 = tau_m.*((b(j,:).*dvdy).*t_m + (yres).*b(j,:).*dbdy(k,:));
%             s_va2 = tau_m.*((gamma-1)/2*b(j,:).*dady).*t_m;
%             
%             % integrate the residuals function
%             vu1_val(count) = sum(weight.*(g_vu1));         vu2_val(count) = sum(weight.*(s_vu2 + g_vu2)); 
%             vv1_val(count) = sum(weight.*(s_vv1 + g_vv1)); vv2_val(count) = sum(weight.*(s_vv2 + g_vv2));
%             va1_val(count) = sum(weight.*(s_va1 + g_va1)); va2_val(count) = sum(weight.*(s_va2 + g_va2));
%             
%             %% calculate the p block (p test functions)
%             % get the stabilization contribution (Part 1/continuity residual
%             g_au1 = tau_c.*((gamma-1)/2*a.*dbdx(j,:)).*c_mc;
%             g_av1 = tau_c.*((gamma-1)/2*a.*dbdy(j,:)).*c_mc;
%             g_aa1 = tau_c.*(u.*dbdx(j,:) + v.*dbdy(j,:)).*c_mc;
%             
%             % get the grad div stabilization (Part 2/continuity residual)
%             g_au2 = tau_c.*((b(j,:).*dadx).*c_mc + (b(j,:).*(cres)).*dbdx(k,:));
%             g_av2 = tau_c.*((b(j,:).*dady).*c_mc + (b(j,:).*(cres)).*dbdy(k,:));
%             g_aa2 = tau_c.*((gamma-1)/2*b(j,:).*(dudx + dvdy)).*c_mc;
%             
%             % get the stabilization contribution (Part 1/ Oseen Part)
%             s_au1 = tau_m.*(-mue*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*dbdx(k,:); % 
%             s_av1 = tau_m.*(-mue*(d2bdx2(j,:) + d2bdy2(j,:)) + (u.*dbdx(j,:) + v.*dbdy(j,:))).*dbdy(k,:); % 
%             s_aa1 = tau_m.*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:));
%             
%             % get the stabilization contribution (Part 2)
%             s_au2 = tau_m.*((b(j,:).*dudx.*dbdx(k,:) + b(j,:).*dvdx.*dbdy(k,:)));
%             s_av2 = tau_m.*((b(j,:).*dudy.*dbdx(k,:) + b(j,:).*dvdy.*dbdy(k,:)));
%             
%             % integrate the residuals function
%             au1_val(count) = sum(weight.*(s_au1 + g_au1)); au2_val(count) = sum(weight.*(s_au2 + g_au2));
%             av1_val(count) = sum(weight.*(s_av1 + g_av1)); av2_val(count) = sum(weight.*(s_av2 + g_av2));
%             aa1_val(count) = sum(weight.*(s_aa1 + g_aa1)); aa2_val(count) = sum(weight.*(g_aa2));
%             
%             %% do other stuff
%             % get the global indices for the assembly
%             j_idx(count) = curr_elem(j); k_idx(count) = curr_elem(k);
%             
%             % count one up
%             count = count + 1;
%}