function [N1, N2, Nres] = get_OP_isaNavierStokes_CTA(self, curr_cond, ...
    red_elems, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the incompressible navier stokes

% check the systems variables (u, v, a)
if sum(ismember(self.sys_variables, {'u', 'v', 'a'})) ~= 3
    error('this functions works for the field variables u, v, and a')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; gamma = 1.4;
n = size(self.mesh.nodes,1); [n_elem, n_base] = size(self.mesh.elems);
n_red_elem = size(red_elems, 1);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_red_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu1_val = temp; uu2_val = temp; uv2_val = temp; ua1_val = temp; ua2_val = temp;
vv1_val = temp; vv2_val = temp; vu2_val = temp; va1_val = temp; va2_val = temp;
au1_val = temp; av1_val = temp; aa1_val = temp;
au2_val = temp; av2_val = temp; aa2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], 3);

% get the starting points of each block and sort them into vectors
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n; uvec = curr_cond(:,uid);
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n; vvec = curr_cond(:,vid);
aid = find(ismember(self.sys_variables,'a')); as = (aid-1) * n; avec = curr_cond(:,aid);

% get the turbulence field for the calculation
if isfield(opts, 'mtvec'), mtvec = opts.mtvec;
else, mtvec = zeros(size(uvec)); end

% loop over all the elements
count = 1;
for i=1:n_red_elem
    
    %% preallocate the mappings of the current element
    % get the current element and its mappings
    curr_elem = red_elems(i,:);
    
    % get shortcuts to the mapping functions
    weight = self.weights(i,:);
    dbdx = self.dbdx(:,:,i); dbdy = self.dbdy(:,:,i);
    uval = uvec(curr_elem); vval = vvec(curr_elem); aval = avec(curr_elem);
    mtval = mtvec(curr_elem);
    
    % get the current velocitys at the quad points
    u = sum(b.*uval,1); dudx = sum(dbdx.*uval,1); dudy = sum(dbdy.*uval,1);
    v = sum(b.*vval,1); dvdx = sum(dbdx.*vval,1); dvdy = sum(dbdy.*vval,1);
    a = sum(b.*aval,1); dadx = sum(dbdx.*aval,1); dady = sum(dbdy.*aval,1);
    mt = sum(b.*mtval,1); % get the turbulent viscosity at quad points
    
    % loop over the base functions twice to use them as test and trial
    for k=1:n_base % loop over the test functions
        for j=1:n_base % loop over the shape functions
            
            %% calculate the u block
            % get the galerkin contrubution
            r_uu1 = (mue + mt).*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                (u.*dbdx(j,:) + v.*dbdy(j,:)).*b(k,:); % convective part
            r_uu2 = (b(j,:).*dudx).*b(k,:); r_uv2 = (b(j,:).*dudy).*b(k,:); % convective part from linearization
            r_ua1 = (2/(gamma-1)*a.*dbdx(j,:)).*b(k,:); % "pressure" contribution
            r_ua2 = (2/(gamma-1)*b(j,:).*dadx).*b(k,:); % "pressure" contribution
            
            % integrate the residuals function
            uu1_val(count) = sum(weight.*(r_uu1));
            uu2_val(count) = sum(weight.*(r_uu2)); uv2_val(count) = sum(weight.*(r_uv2));
            ua1_val(count) = sum(weight.*(r_ua1)); ua2_val(count) = sum(weight.*(r_ua2));
            
            %% calculate the v block
            % get the galerkin contrubution
            r_vu2 = (b(j,:).*dvdx).*b(k,:); r_vv2 = (b(j,:).*dvdy).*b(k,:); % convective part from linearization
            r_vv1 = (mue + mt).*(dbdx(j,:).*dbdx(k,:) + dbdy(j,:).*dbdy(k,:)) +... % diffusive part
                (u.*dbdx(j,:) + v.*dbdy(j,:)).*b(k,:); %convective  part
            r_va1 = (2/(gamma-1)*a.*dbdy(j,:)).*b(k,:); % "pressure" contribution
            r_va2 = (2/(gamma-1)*b(j,:).*dady).*b(k,:);
            
            % integrate the residuals function
            vu2_val(count) = sum(weight.*(r_vu2)); vv2_val(count) = sum(weight.*(r_vv2));
            vv1_val(count) = sum(weight.*(r_vv1));
            va1_val(count) = sum(weight.*(r_va1)); va2_val(count) = sum(weight.*(r_va2));
            
            %% calculate the p block
            % get the galerkin contrubution
            r_au1 = (gamma-1)/2*(a.*dbdx(j,:)).*b(k,:); r_au2 = (b(j,:).*dadx).*b(k,:);
            r_av1 = (gamma-1)/2*(a.*dbdy(j,:)).*b(k,:); r_av2 = (b(j,:).*dady).*b(k,:);
            r_aa1 = (u.*dbdx(j,:) + v.*dbdy(j,:)).*b(k,:); 
            r_aa2 = (gamma-1)/2*(b(j,:).*(dudx + dvdy)).*b(k,:);
            
            % integrate the residuals function
            au1_val(count) = sum(weight.*(r_au1)); au2_val(count) = sum(weight.*(r_au2));
            av1_val(count) = sum(weight.*(r_av1)); av2_val(count) = sum(weight.*(r_av2));
            aa1_val(count) = sum(weight.*(r_aa1)); aa2_val(count) = sum(weight.*(r_aa2));
            
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
