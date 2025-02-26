function [opts, curr_cond] = get_Turbulence_field_SA(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the incompressible navier stokes

% check the systems variables (u, v, p)
if sum(ismember(self.sys_variables, {'u', 'v', 'nt'})) ~= 3
    error('this functions works for the field variables u, v, and nt')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; rho = opts.rho;
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); 
j_idx = temp; k_idx = temp;
nn1_val = temp; nn2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);
curr_cond = reshape(curr_cond, [], numel(self.sys_variables));

% get the u velocities and its derivatives at the int points
uvec = curr_cond(:,ismember(self.sys_variables,'u'));
u = self.mesh.get_fun_val(uvec, 'int'); 
[~, dudy] = self.mesh.get_1st_deriv(uvec, 'int');

% get the v velocities and its derivatives at the int points
vvec = curr_cond(:,ismember(self.sys_variables,'v'));
v = self.mesh.get_fun_val(vvec, 'int'); 
[dvdx, ~] = self.mesh.get_1st_deriv(vvec, 'int');

% get the turbulent viscosity at the int points
ntvec = curr_cond(:,ismember(self.sys_variables,'nt'));
nt = self.mesh.get_fun_val(ntvec, 'int'); 
[dntdx, dntdy] = self.mesh.get_1st_deriv(ntvec, 'int');
d = self.mesh.get_fun_val(self.wall_dist, 'int');

%% calculate some stuff for the turbulence model
% specify the constants for the turbulence model
cb1 = 0.1355; cb2 = 0.622; cv1 = 7.1; sig = 2/3; kappa = 0.41;
cw1 = cb1/kappa^2 + (1+cb2)/sig; cw2 = 0.3; cw3 = 2;
cv2  = 0.7; cv3 = 0.9; cn1 = 16;

% calculate some values for the SA turbulence model
X  = nt./(mue/rho);
fv1 = X.^3 ./ (X.^3 + cv1^3);
fv2 = 1 - X ./ (1 + X.*fv1);

% get the modified version of the strain rate to prevent negative
% values of Se (Pepijn S. 16)
S = sqrt((dudy - dvdx).^2); % strain rate within the flowfield
Sm = nt./(kappa^2*d.^2).*fv2;
Smm = (S.*(cv2^2*S + cv3*Sm))./...
    (S.*(cv3 - 2*cv2) - Sm);
Se = S + Sm; sel = Sm < -cv2*S; 
Se(sel) = S(sel) + Smm(sel);

% calculate the other functions
r = nt./(kappa^2*d.^2.*Se); r(r>10) = 10;
g = r + cw2*(r.^6 - r);
fw = g.*((1 + cw3^6)./(g.^6 + cw3^6)).^(1/6);

% negative SA-model
fn = (cn1 + X.^3)./(cn1 - X.^3);

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
            %% calculate the nue_t block
            % get the galerkin contribution
            r_nn1 = (u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:)).*b(k,:) + ...  % the convective part
                - cb1.*Se(i,:).*b(j,:).*b(k,:) + ...
                + cw1.*fw(i,:).*nt(i,:)./d(i,:).^2.*b(j,:).*b(k,:) + ...
                + 1/sig*(mue/rho + nt(i,:)).*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) + ...
                - 1/sig*cb2*(dntdx(i,:).*dbdx(j,:) + dntdy(i,:).*dbdy(j,:)).*b(k,:);
            r_nn2 = + cw1.*fw(i,:).*nt(i,:)./d(i,:).^2.*b(j,:).*b(k,:);
            
            % integrate the residual function
            nn1_val(count) = r_nn1*weight; nn2_val(count) = r_nn2*weight;
            
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
N1 = sparse(k_idx, j_idx, nn1_val, n, n);
N2 = sparse(k_idx, j_idx, nn2_val, n, n);

% get the residuals at the nodes and the linearized matrix
Nres = N1*ntvec; N = N1 + N2;

% get the dirichlet bounds of the turbulence field
dir_sel = reshape(self.dir_bounds(:,1), [], numel(self.sys_variables));
i_DoF = ~dir_sel(:, ismember(self.sys_variables, 'nt'));

% solve the system
ntvec(i_DoF) = ntvec(i_DoF) - N(i_DoF, i_DoF)\Nres(i_DoF);
curr_cond(:, ismember(self.sys_variables,'nt')) = ntvec;
curr_cond = reshape(curr_cond, [], 1);

% get the field of turbulence viscosity
Xvec  = ntvec./(mue/rho);
fv1vec = Xvec.^3./(Xvec.^3 + cv1^3);
opts.mtvec = ntvec.*fv1vec * rho;
end
