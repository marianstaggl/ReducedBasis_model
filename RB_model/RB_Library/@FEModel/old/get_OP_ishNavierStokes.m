function [N1, N2, Nres] = get_OP_isNavierStokes(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the incompressible navier stokes

% check the systems variables (u, v, p)
if sum(ismember(self.sys_variables, {'u', 'v', 'h'})) ~= 3
    error('this functions works for the field variables u, v, and p')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; gamma = 1.4;
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; uu2_val = temp; uv2_val = temp; uh_val = temp;
vv1_val = temp; vv2_val = temp; vu2_val = temp; vh_val = temp;
hu1_val = temp; hv1_val = temp; hh1_val = temp;
hu2_val = temp; hv2_val = temp; hh2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], 3);

% get the starting points of each block and sort them into vectors
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);
hid = find(ismember(self.sys_variables,'h')); hs = (hid-1) * n; hvec = curr_cond(:,hid);

% loop over all the elements
count = 1;
for i=1:n_elem
    
    %% preallocate the mappings of the current element
    % get the current element and its mappings
    curr_elem = self.mesh.elems(i,:);
    
    % get shortcuts to the mapping functions
    weight = self.weights(i,:);
    dbdx = self.dbdx(:,:,i); dbdy = self.dbdy(:,:,i);
    uval = uvec(curr_elem); vval = vvec(curr_elem); hval = hvec(curr_elem);
    
    % get the current velocitys at the quad points
    u = sum(b.*uval,1); dudx = sum(dbdx.*uval,1); dudy = sum(dbdy.*uval,1);
    v = sum(b.*vval,1); dvdx = sum(dbdx.*vval,1); dvdy = sum(dbdy.*vval,1);
    h = sum(b.*hval,1); dhdx = sum(dbdx.*hval,1); dhdy = sum(dbdy.*hval,1);
    
    % loop over the base functions twice to use them as test and trial
    for k=1:n_base % loop over the test functions
        for j=1:n_base % loop over the shape functions
            
            %% calculate the u block
            % get the galerkin contrubution
            r_uu1 = mue*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                (u.*dbdx(j,:) + v.*dbdy(j,:)).*b(k,:); % convective part
            r_uu2 = (b(j,:).*dudx).*b(k,:); r_uv2 = (b(j,:).*dudy).*b(k,:); % convective part from linearization
            r_uh = dbdx(j,:).*b(k,:); % pressure contribution
            
            % integrate the residuals function
            uu1_val(count) = sum(weight.*(r_uu1));
            uu2_val(count) = sum(weight.*(r_uu2)); uv2_val(count) = sum(weight.*(r_uv2));
            uh_val(count) = sum(weight.*(r_uh));
            
            %% calculate the v block
            % get the galerkin contrubution
            r_vu2 = (b(j,:).*dvdx).*b(k,:); r_vv2 = (b(j,:).*dvdy).*b(k,:); % convective part from linearization
            r_vv1 = mue*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                (u.*dbdx(j,:) + v.*dbdy(j,:)).*b(k,:); %convective  part
            r_vh = dbdy(j,:).*b(k,:); % pressure contribution
            
            % integrate the residuals function
            vu2_val(count) = sum(weight.*(r_vu2)); vv2_val(count) = sum(weight.*(r_vv2));
            vv1_val(count) = sum(weight.*(r_vv1));
            vh_val(count) = sum(weight.*(r_vh));
            
            %% calculate the p block
            % get the galerkin contrubution
            r_hu1 = (gamma-1)*(h.*dbdx(j,:)).*b(k,:); r_hu2 = (b(j,:).*dhdx).*b(k,:);
            r_hv1 = (gamma-1)*(h.*dbdy(j,:)).*b(k,:); r_hv2 = (b(j,:).*dhdy).*b(k,:);
            r_hh1 = (u.*dbdx(j,:) + v.*dbdy(j,:)).*b(k,:); r_hh2 = (gamma-1)*(b(j,:).*(dudx + dvdy)).*b(k,:);
            
            % integrate the residuals function
            hu1_val(count) = sum(weight.*(r_hu1)); hu2_val(count) = sum(weight.*(r_hu2));
            hv1_val(count) = sum(weight.*(r_hv1)); hv2_val(count) = sum(weight.*(r_hv2));
            hh1_val(count) = sum(weight.*(r_hh1)); hh2_val(count) = sum(weight.*(r_hh2));
            
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
v1_glob = [uu1_val,    uh_val,     vv1_val,    vh_val,     hu1_val,    hv1_val,    hh1_val   ];
k1_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + hs, k_idx + hs, k_idx + hs];
j1_glob = [j_idx + us, j_idx + hs, j_idx + vs, j_idx + hs, j_idx + us, j_idx + vs, j_idx + hs];

% reshape the vetors for the sparse matrix assembly
v2_glob = [uu2_val,    uv2_val,    vu2_val,    vv2_val,    hu2_val,    hv2_val,    hh2_val   ];
k2_glob = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + hs, k_idx + hs, k_idx + hs];
j2_glob = [j_idx + us, j_idx + vs, j_idx + us, j_idx + vs, j_idx + us, j_idx + vs, j_idx + vs];

% get a sparse global matrix
N1 = sparse(k1_glob, j1_glob, v1_glob, numel(curr_cond), numel(curr_cond));
N2 = sparse(k2_glob, j2_glob, v2_glob, numel(curr_cond), numel(curr_cond));

% get the residuals at the nodes and the linearized matrix
Nres = N1*reshape(curr_cond,[],1);

end
