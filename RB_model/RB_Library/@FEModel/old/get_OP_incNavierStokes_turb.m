function [N1, N2, Nres] = get_OP_incNavierStokes_turb(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the incompressible navier stokes

% check the systems variables (u, v, p, nt)
if sum(ismember(self.sys_variables, {'u', 'v', 'p', 'nt'})) ~= 4
    error('this functions works for the field variables u, v, p and nt')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; rho = opts.rho;
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; uu2_val = temp; uv2_val = temp; up_val = temp;
vv1_val = temp; vv2_val = temp; vu2_val = temp; vp_val = temp;
pu1_val = temp; pv1_val = temp; nn1_val = temp; nn2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], 4);

% get the starting points of each block and sort them into vectors
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);
pid = find(ismember(self.sys_variables,'p')); ps = (pid-1) * n; pvec = curr_cond(:,pid);
ntid = find(ismember(self.sys_variables,'nt')); nts = (ntid-1) * n; ntvec = curr_cond(:,ntid);
yvec = self.wall_dist;

% specify the constants for the turbulence model
c.sig = 2/3; c.cb1 = 0.1355; c.cb2 = 0.622; c.kappa = 0.41;
c.cw1 = c.cb1/c.kappa^2 + (1+c.cb2)/c.sig; c.cw2 = 0.3; c.cw3 = 2;
c.cv1 = 7.1; c.ct1 = 1; c.ct2 = 2; c.ct3 = 1.2; c.ct4 = 0.5;
c.rho = rho; c.mue = mue;

% loop over all the elements
count = 1;
for i=1:n_elem
    
    %% preallocate the mappings of the current element
    % get the current element and its mappings
    curr_elem = self.mesh.elems(i,:);
    
    % get shortcuts to the mapping functions
    weight = self.weights(i,:);  yval = yvec(curr_elem);
    dbdx = self.dbdx(:,:,i); dbdy = self.dbdy(:,:,i);
    uval = uvec(curr_elem); vval = vvec(curr_elem); 
    pval = pvec(curr_elem); ntval = ntvec(curr_elem);
    
    % get the current velocitys at the quad points
    u = sum(b.*uval,1); dudx = sum(dbdx.*uval,1); dudy = sum(dbdy.*uval,1);
    v = sum(b.*vval,1); dvdx = sum(dbdx.*vval,1); dvdy = sum(dbdy.*vval,1);
    nt = sum(b.*ntval,1); dntdx = sum(dbdx.*ntval,1); dntdy = sum(dbdy.*ntval,1);
    y = sum(b.*yval,1); % get the wall distance at the quad points;
    
    % get the turbulent values at the quad points
    [mt, C1, C2, C3, C4] = get_turb_vals(dudy, dvdx, nt, y, c);
    
    % loop over the base functions twice to use them as test and trial
    for k=1:n_base % loop over the test functions
        for j=1:n_base % loop over the shape functions
            %% calculate the u block
            % get the galerkin contrubution
            r_uu1 = (mue + mt) .*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                rho*(u.*dbdx(j,:) + v.*dbdy(j,:)).*b(k,:); % convective part
            r_uu2 = rho*(b(j,:).*dudx).*b(k,:); r_uv2 = rho*(b(j,:).*dudy).*b(k,:); % convective part from linearization
            r_up = dbdx(j,:).*b(k,:); % pressure contribution (integrated by parts) old --> - b(j,:).*dbdx(k,:); % pressure contribution
            
            % integrate the residuals function
            uu1_val(count) = sum(weight.*(r_uu1));
            uu2_val(count) = sum(weight.*(r_uu2)); uv2_val(count) = sum(weight.*(r_uv2));
            up_val(count) = sum(weight.*(r_up));
            
            %% calculate the v block
            % get the galerkin contrubution
            r_vu2 = rho*(b(j,:).*dvdx).*b(k,:); r_vv2 = rho*(b(j,:).*dvdy).*b(k,:); % convective part from linearization
            r_vv1 = (mue + mt) .*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                rho*(u.*dbdx(j,:) + v.*dbdy(j,:)).*b(k,:); %convective  part
            r_vp = dbdy(j,:).*b(k,:); % pressure contribution (integrated by parts) old --> - b(j,:).*dbdy(k,:);
            
            % integrate the residuals function
            vu2_val(count) = sum(weight.*(r_vu2)); vv2_val(count) = sum(weight.*(r_vv2));
            vv1_val(count) = sum(weight.*(r_vv1));
            vp_val(count) = sum(weight.*(r_vp));
            
            %% calculate the p block
            % get the galerkin contrubution
            r_pu = dbdx(j,:).*b(k,:);
            r_pv = dbdy(j,:).*b(k,:);
            
            % integrate the residuals function
            pu1_val(count) = sum(weight.*(r_pu));
            pv1_val(count) = sum(weight.*(r_pv));
            
            %% calculate the nue_t block
            % get the galerkin contribution
            r_nn1 = (u.*dbdx(j,:) + v.*dbdy(j,:) - C1.*b(j,:) +...
                C2.*nt./(y.^2).*b(j,:) - C4.*(dntdx.*dbdx(j,:) +...
                dntdy.*dbdy(j,:))).*b(k,:) + C3.*(dbdx(j,:).*dbdx(k,:) +...
                dbdx(j,:).*dbdy(k,:));
            r_nn2 = ( C2.*nt./(y.^2).*b(j,:) - C4.*(dntdx.*dbdx(j,:) +...
                dntdy.*dbdy(j,:))).*b(k,:);
            
            % integrate the residual function
            nn1_val(count) = sum(weight.*(r_nn1)); nn2_val(count) = sum(weight.*(r_nn2));
            
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
v1_glob = [uu1_val,    up_val,     vv1_val,    vp_val,     pu1_val,     pv1_val,   nn1_val    ];
k1_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + ps, k_idx + ps, k_idx + nts];
j1_glob = [j_idx + us, j_idx + ps, j_idx + vs, j_idx + ps, j_idx + us, j_idx + vs, j_idx + nts];

% reshape the vetors for the sparse matrix assembly
v2_glob = [uu2_val,    uv2_val,    vu2_val,    vv2_val];
k2_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs];
j2_glob = [j_idx + us, j_idx + vs, j_idx + us, j_idx + vs];

% reshape the vectors for the sparse matrix assembly
v3_glob = nn2_val;
k3_glob = k_idx + nts;
j3_glob = j_idx + nts;

% get a sparse global matrix
N1 = sparse(k1_glob, j1_glob, v1_glob, numel(curr_cond), numel(curr_cond));
N2 = sparse(k2_glob, j2_glob, v2_glob, numel(curr_cond), numel(curr_cond));
N3 = sparse(k3_glob, j3_glob, v3_glob, numel(curr_cond), numel(curr_cond));

% get the residuals at the nodes and the linearized matrix
Nres = N1*reshape(curr_cond,[],1); N1 = N1 + N3;
end

function [mt, C1, C2, C3, C4] = get_turb_vals(dudy, dvdx, nt, y, c)
% calculate the "constants" for the SA iter
X  = nt./(c.mue./c.rho); 
fv1 = X.^3./(X.^3 + c.cv1^3);
fv2 = 1 - X./(1 + X.*fv1);
% ft2 = c.ct3*exp(-c.ct4*X.^2);
omega_n = sqrt((dudy - dvdx).^2);
st = omega_n + nt./(c.kappa^2*y.^2).*fv2;
r = nt./(c.kappa^2*st.*y.^2); r(r>10) = 10;
g = r + c.cw2*(r.^6 - r);
fw = g.*((1+c.cw3^6)./(g.^6 + c.cw3^6)).^(1/6);

% group them together to new "constants"
C1 = c.cb1.*st; %c.cb1.*(1-ft2).*st;
C2 = c.cw1*fw; % (c.cw1*fw - c.cb1/c.kappa*ft2);
C3 = 1/c.sig*(c.mue./c.rho + nt);
C4 = c.cb2/c.sig;

% calculate the turbulent viscosity
mt = c.rho.*nt.*fv1;
end
