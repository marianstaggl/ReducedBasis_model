function [N, Nres] = get_OP_StokesFlow(self, curr_cond, opts)
%GET_OP_ISNAVIERSTOKES get the operator for the incompressible navier stokes

% check the systems variables (u, v, p)
if sum(ismember(self.sys_variables, {'u', 'v', 'p'})) ~= 3
    error('this functions works for the field variables u, v, and p')
end

% get the basefunction at the integration points
n = size(self.mesh.nodes,1);
int_points = self.mesh.ref_elem.int_points;
b = self.mesh.ref_elem.get_b(int_points);

% get the locations of the different variables
uid = find(ismember(self.sys_variables,'u')); us = (uid-1) * n;
vid = find(ismember(self.sys_variables,'v')); vs = (vid-1) * n;
pid = find(ismember(self.sys_variables,'p')); ps = (pid-1) * n;

% build the systems matrix
% [k_idx, j_idx, s_val] = build_StokesSys(self.mesh.elems, b, self.dbdx, self.dbdy, self.weights, opts.mue, us, vs, ps);
[k_idx, j_idx, s_val] = self.build_StokesSys(self.mesh.elems, b, self.dbdx, self.dbdy, self.weights, opts.mue, us, vs, ps);


% get the residuals at the nodes and the linearized matrix
N = sparse(k_idx, j_idx, s_val, numel(curr_cond), numel(curr_cond));
Nres = N*reshape(curr_cond,[],1);
end

function [k_idx, j_idx, s_val] = build_StokesSys(elems, b, dbdx, dbdy, weights, mue, us, vs, ps)

% get the number of elements and base functions
[n_elem, n_base] = size(elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1, (n_elem*n_base*n_base)); j_idx = temp; k_idx = temp;
uu_val = temp; up_val = temp; vv_val = temp; vp_val = temp;
pu_val = temp; pv_val = temp; pp_val = temp;


% loop over all the elements
count = 1;
for i=1:n_elem
    
    %% preallocate the mappings of the current element
    % get the current element and its mappings
    curr_elem = elems(i,:);
    
    % get shortcuts to the mapping functions
    c_weight = weights(i,:)'; dbdx_c = dbdx(:,:,i); dbdy_c = dbdy(:,:,i);
    
    % loop over the base functions twice to use them as test and trial
    for k=1:n_base % loop over the test functions
        for j=1:n_base % loop over the shape functions
            
            %% calculate the u block
            % get the galerkin contrubution
            r_uu = mue*(dbdx_c(j,:).*dbdx_c(k,:) + dbdy_c(j,:).*dbdy_c(k,:));
            r_up = dbdx_c(j,:).*b(k,:); % pressure contribution (integrated by parts) old --> - b(j,:).*dbdx_c(k,:);
            
            % integrate the residuals function
            uu_val(count) = r_uu*c_weight;
            up_val(count) = r_up*c_weight;
            
            %% calculate the v block
            % get the galerkin contrubution
            r_vv = mue*(dbdx_c(j,:).*dbdx_c(k,:) + dbdy_c(j,:).*dbdy_c(k,:));
            r_vp = dbdy_c(j,:).*b(k,:); % pressure contribution (integrated by parts) old --> - b(j,:).*dbdy_c(k,:);
            
            % integrate the residuals function
            vv_val(count) = r_vv*c_weight;
            vp_val(count) = r_vp*c_weight;
            
            %% calculate the p block
            % get the galerkin contrubution
            r_pu = dbdx_c(j,:).*b(k,:);
            r_pv = dbdy_c(j,:).*b(k,:);
            r_pp = (dbdx_c(j,:).*dbdx_c(k,:) + dbdy_c(j,:).*dbdy_c(k,:));
            
            % integrate the residuals function
            pu_val(count) = r_pu*c_weight;
            pv_val(count) = r_pv*c_weight;
            pp_val(count) = r_pp*c_weight;
            
            %% do other stuff
            % get the global indices for the assembly
            j_idx(count) = curr_elem(j); k_idx(count) = curr_elem(k);
            
            % count one up
            count = count + 1;
        end
    end
end

% reshape the vectors for the sparse matrix assembly
s_val = [uu_val,     up_val,     vv_val,     vp_val,     pu_val,     pv_val    , pp_val];
k_idx = [k_idx + us, k_idx + us, k_idx + vs, k_idx + vs, k_idx + ps, k_idx + ps, k_idx + ps];
j_idx = [j_idx + us, j_idx + ps, j_idx + vs, j_idx + ps, j_idx + us, j_idx + vs, j_idx + ps];
end
