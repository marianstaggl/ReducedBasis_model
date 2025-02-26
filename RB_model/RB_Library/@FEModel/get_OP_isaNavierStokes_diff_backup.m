function [N1_diff, R_diff] = get_OP_isaNavierStokes_diff_backup(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the incompressible navier stokes

% check the systems variables (u, v, a)
if sum(ismember(self.sys_variables, {'u', 'v', 'a'})) ~= 3
    error('this functions works for the field variables u, v, and a')
end

%% preallocate various variables
% get shortcuts to the physical properties and cell counts
mue = opts.mue; n = size(self.mesh.nodes,1); 
[n_elem, n_base] = size(self.mesh.elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); 
j_idx = temp; k_idx = temp;
uu1_val = temp; vv1_val = temp;
uu2_val = temp; vv2_val = temp;

% get the basefunction at the integration points
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the current conditions and reshape them into a new order
curr_cond = reshape(curr_cond, [], 3);

% get the starting points of each block and sort them into vectors
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n;
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n;

% get the turbulence field for the calculation
if isfield(opts, 'mtvec'), mtvec = opts.mtvec;
else, mtvec = zeros(size(curr_cond(:,uid))); end

% loop over all the elements
count = 1;
for i=1:n_elem
    
    %% preallocate the mappings of the current element
    % get the current element and its mappings
    curr_elem = self.mesh.elems(i,:);
    
    % get shortcuts to the mapping functions
    weight = self.weights(i,:);
    dbdx = self.dbdx(:,:,i); 
    dbdy = self.dbdy(:,:,i);
    mtval = mtvec(curr_elem);
    
    % get the current velocitys at the quad points
    mt = sum(b.*mtval,1); 
    dmtdx = sum(dbdx.*mtval,1); 
    dmtdy = sum(dbdy.*mtval,1); % get the turbulent viscosity at quad points
    
    % loop over the base functions twice to use them as test and trial
    for k=1:n_base % loop over the test functions
        for j=1:n_base % loop over the shape functions
            
            %% calculate the u block
            % get the galerkin contrubution
            r_uu1 = (mue + mt).*(dbdx(j,:).*dbdx(k,:) +...
                dbdy(j,:).*dbdy(k,:)); % diffusive part
            
            % integrate the residuals function
            uu1_val(count) = sum(weight.*(r_uu1));
            
            %% calculate the v block
            % get the galerkin contrubution
            r_vv1 = (mue + mt).*(dbdx(j,:).*dbdx(k,:) +...
                dbdy(j,:).*dbdy(k,:)); % diffusive part
            
            % integrate the residuals function
            vv1_val(count) = sum(weight.*(r_vv1));
            
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
v1_glob = [uu1_val,    vv1_val   ];
k1_glob = [k_idx + us, k_idx + vs];
j1_glob = [j_idx + us, j_idx + vs];

% get a sparse global matrix
N1_diff = sparse(k1_glob, j1_glob, v1_glob, numel(curr_cond), numel(curr_cond));
R_diff = N1_diff*reshape(curr_cond,[],1);

end
