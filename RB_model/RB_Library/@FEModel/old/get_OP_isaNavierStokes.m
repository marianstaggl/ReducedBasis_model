function [N1, N2, Nres] = get_OP_isaNavierStokes(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the incompressible navier stokes

% check the systems variables (u, v, a)
if sum(ismember(self.sys_variables, {'u', 'v', 'a'})) ~= 3
    error('this functions works for the field variables u, v, and a')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; gamma = 1.4;
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; uu2_val = temp; uv2_val = temp; ua1_val = temp; ua2_val = temp;
vv1_val = temp; vv2_val = temp; vu2_val = temp; va1_val = temp; va2_val = temp;
au1_val = temp; av1_val = temp; aa1_val = temp;
au2_val = temp; av2_val = temp; aa2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], numel(self.sys_variables));

% get the starting points of each block and sort them into vectors
% get the u velocities and its derivatives at the int points
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
u = self.mesh.get_fun_val(uvec, 'int'); [dudx, dudy] = self.mesh.get_1st_deriv(uvec, 'int');

% get the v velocities and its derivatives at the int points
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);
v = self.mesh.get_fun_val(vvec, 'int'); [dvdx, dvdy] = self.mesh.get_1st_deriv(vvec, 'int');

% get the pressure its derivatives at the int points
aid = find(ismember(self.sys_variables,'a')); as = (aid-1) * n; avec = curr_cond(:,aid);
a = self.mesh.get_fun_val(avec, 'int'); [dadx, dady] = self.mesh.get_1st_deriv(avec, 'int');

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
            r_uu1 = mue * (dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                (u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:)).*b(k,:); % convective part
            r_uu2 = (b(j,:).*dudx(i,:)).*b(k,:); r_uv2 = (b(j,:).*dudy(i,:)).*b(k,:); % convective part from linearization
            r_ua1 = (2/(gamma-1)*a(i,:).*dbdx(j,:)).*b(k,:); % with integration by parts: -(1/(gamma-1)*b(k,:).*a(i,:).*dbdx(k,:));
            r_ua2 = (2/(gamma-1)*b(j,:).*dadx(i,:)).*b(k,:); % with integration by parts: -(1/(gamma-1)*a(i,:).*b(k,:).*dbdx(k,:));
            
            % integrate the residuals function
            uu1_val(count) = r_uu1*weight;
            uu2_val(count) = r_uu2*weight; uv2_val(count) = r_uv2*weight;
            ua1_val(count) = r_ua1*weight; ua2_val(count) = r_ua2*weight;
            
            %% calculate the v block
            % get the galerkin contrubution
            r_vu2 = (b(j,:).*dvdx(i,:)).*b(k,:); r_vv2 = (b(j,:).*dvdy(i,:)).*b(k,:); % convective part from linearization
            r_vv1 = mue * (dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                (u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:)).*b(k,:); %convective  part
            r_va1 = (2/(gamma-1)*a(i,:).*dbdy(j,:)).*b(k,:); % with integration by parts: -(1/(gamma-1)*b(k,:).*a(i,:).*dbdy(k,:));
            r_va2 = (2/(gamma-1)*b(j,:).*dady(i,:)).*b(k,:); % with integration by parts: -(1/(gamma-1)*a(i,:).*b(k,:).*dbdy(k,:));
            
            % integrate the residuals function
            vu2_val(count) = r_vu2*weight; vv2_val(count) = r_vv2*weight;
            vv1_val(count) = r_vv1*weight;
            va1_val(count) = r_va1*weight; va2_val(count) = r_va2*weight;
            
            %% calculate the p block
            % get the galerkin contrubution
            r_au1 = (gamma-1)/2*(a(i,:).*dbdx(j,:)).*b(k,:); r_au2 = (b(j,:).*dadx(i,:)).*b(k,:);
            r_av1 = (gamma-1)/2*(a(i,:).*dbdy(j,:)).*b(k,:); r_av2 = (b(j,:).*dady(i,:)).*b(k,:);
            r_aa1 = (u(i,:).*dbdx(j,:) + v(i,:).*dbdy(j,:)).*b(k,:); 
            r_aa2 = (gamma-1)/2*(b(j,:).*(dudx(i,:) + dvdy(i,:))).*b(k,:);
            
            % integrate the residuals function
            au1_val(count) = r_au1*weight; au2_val(count) = r_au2*weight;
            av1_val(count) = r_av1*weight; av2_val(count) = r_av2*weight;
            aa1_val(count) = r_aa1*weight; aa2_val(count) = r_aa2*weight;
            
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

% reshape the vetors for the sparse matrix assembly
v2_glob = [uu2_val,    uv2_val,    ua2_val,    vu2_val,    vv2_val,    va2_val,    au2_val,    av2_val,    aa2_val   ];
k2_glob = [k_idx + us, k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + vs, k_idx + as, k_idx + as, k_idx + as];
j2_glob = [j_idx + us, j_idx + vs, j_idx + as, j_idx + us, j_idx + vs, j_idx + as, j_idx + us, j_idx + vs, j_idx + vs];

% get a sparse global matrix
N1 = sparse(k1_glob, j1_glob, v1_glob, numel(curr_cond), numel(curr_cond));
N2 = sparse(k2_glob, j2_glob, v2_glob, numel(curr_cond), numel(curr_cond));

% get the residuals at the nodes and the linearized matrix
Nres = N1*reshape(curr_cond,[],1);

end
