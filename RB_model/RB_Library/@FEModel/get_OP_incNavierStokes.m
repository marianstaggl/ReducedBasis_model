function [N1, N2, Nres] = get_OP_incNavierStokes(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the incompressible navier stokes

% check the systems variables (u, v, p)
if sum(ismember(self.sys_variables, {'u', 'v', 'p'})) ~= 3
    error('this functions works for the field variables u, v, and p')
end

%% preallocate various variables
% get shortcuts to the properties and cell counts
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; uu2_val = temp; uv2_val = temp; up_val = temp;
vv1_val = temp; vv2_val = temp; vu2_val = temp; vp_val = temp;
pu_val = temp; pv_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], numel(self.sys_variables));

% get the u velocities and its derivatives at the int points
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
u = self.mesh.get_fun_val(uvec, 'int'); [dudx, dudy] = self.mesh.get_1st_deriv(uvec, 'int');

% get the v velocities and its derivatives at the int points
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);
v = self.mesh.get_fun_val(vvec, 'int'); [dvdx, dvdy] = self.mesh.get_1st_deriv(vvec, 'int');

% get the pressure its derivatives at the int points
pid = find(ismember(self.sys_variables,'p')); ps = (pid-1) * n; pvec = curr_cond(:,pid);
p = self.mesh.get_fun_val(pvec, 'int');

% get the turbulence viscosity field and physical props
mue = opts.mue; rho = opts.rho;
if ~isfield(opts, 'mtvec'), opts.mtvec = zeros(1, n); end
mt = self.mesh.get_fun_val(opts.mtvec, 'int');

% loop over all the elements
count = 1;
for i=1:n_elem
    
    %% preallocate the mappings of the current element
    % get the current element and its mappings
    curr_elem = self.mesh.elems(i,:);
    
    % get shortcuts to the mapping functions
    weight = self.weights(i,:)'; dbdx = self.dbdx(:,:,i); dbdy = self.dbdy(:,:,i);
    
    % loop over the base functions twice to use them as test and trial
    for k=1:n_base % loop over the test functions
        for j=1:n_base % loop over the shape functions
            %% calculate the u block
            % get the galerkin contrubution
            r_uu1 = (mue + mt(i,:)).*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                rho*(u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:)).*b(k,:); % convective part
            r_uu2 = rho*(b(j,:).*dudx(i,:)).*b(k,:); r_uv2 = rho*(b(j,:).*dudy(i,:)).*b(k,:); % convective part from linearization
            r_up = - b(j,:).*dbdx(k,:); % pressure contribution (integrated by parts) old --> dbdx(j,:).*b(k,:); 
            
            % integrate the residuals function
            uu1_val(count) = r_uu1*weight;
            uu2_val(count) = r_uu2*weight; uv2_val(count) = r_uv2*weight;
            up_val(count) = r_up*weight;
            
            %% calculate the v block
            % get the galerkin contrubution
            r_vu2 = rho*(b(j,:).*dvdx(i,:)).*b(k,:); r_vv2 = rho*(b(j,:).*dvdy(i,:)).*b(k,:); % convective part from linearization
            r_vv1 = (mue + mt(i,:)).*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                rho*(u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:)).*b(k,:); %convective  part
            r_vp = - b(j,:).*dbdy(k,:); % pressure contribution (integrated by parts) old --> dbdy(j,:).*b(k,:); 
            
            % integrate the residuals function
            vu2_val(count) = r_vu2*weight; vv2_val(count) = r_vv2*weight;
            vv1_val(count) = r_vv1*weight;
            vp_val(count) = r_vp*weight;
            
            %% calculate the p block
            % get the galerkin contrubution
            r_pu = dbdx(j,:).*b(k,:);
            r_pv = dbdy(j,:).*b(k,:);
            
            % integrate the residuals function
            pu_val(count) = r_pu*weight;
            pv_val(count) = r_pv*weight;
            
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
v1_glob = [uu1_val,    up_val,     vv1_val,    vp_val,     pu_val,     pv_val    ];
k1_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + ps, k_idx + ps];
j1_glob = [j_idx + us, j_idx + ps, j_idx + vs, j_idx + ps, j_idx + us, j_idx + vs];

% reshape the vetors for the sparse matrix assembly
v2_glob = [uu2_val,    uv2_val,    vu2_val,    vv2_val   ];
k2_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs];
j2_glob = [j_idx + us, j_idx + vs, j_idx + us, j_idx + vs];

% get a sparse global matrix
N1 = sparse(k1_glob, j1_glob, v1_glob, numel(curr_cond), numel(curr_cond));
N2 = sparse(k2_glob, j2_glob, v2_glob, numel(curr_cond), numel(curr_cond));

% get the residuals at the nodes and the linearized matrix
Nres = N1*reshape(curr_cond,[],1);
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
