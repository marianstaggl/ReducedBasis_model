function [L, Lres] = get_OP_PotentialFlow(self, curr_cond, ~)
%GET_OP_HEATTRANSF get the operator for the heat transfer equation

% check the systems variables (this only works for t)
if ~isequal(self.sys_variables, {'phi'})
    error('this function only works for a temperature field!')
end
           
% build the systems matrix
% [k_idx, j_idx, t_val] = self.build_PotentialSys(self.mesh.elems, self.dbdx, self.dbdy, self.weights);
[k_idx, j_idx, t_val] = build_PotentialSys(self.mesh.elems, self.dbdx, self.dbdy, self.weights);

% get the sparse global matrix and the current residual
L = sparse(k_idx, j_idx, t_val); Lres = L*curr_cond;
end

function [k_idx, j_idx, t_val] = build_PotentialSys(elems, dbdx, dbdy, weights)

% preallocate the matrix and set up the numbers
[n_elem, n_base] = size(elems);

% preallocate the vectors for the sparse matrix assembly
temp = zeros(1,n_elem*n_base*n_base);
t_val = temp; j_idx = temp; k_idx = temp;

% loop over all the elements
count = 1;
for i=1:n_elem
    
    % get the current element and its mappings
    curr_elem = elems(i,:); 
    
    % get shortcuts to the mapped derivatives
    c_weight = weights(i,:)'; c_dbdx = dbdx(:,:,i); c_dbdy = dbdy(:,:,i);

    % loop over all of the elements nodes twice
    for k=1:n_base
        for j=1:n_base
            
            % get the real indices and the mapping for the derivatives
            j_idx(count) = curr_elem(j); k_idx(count) = curr_elem(k);
            
            % get the residuals function
            r = (c_dbdx(j,:).*c_dbdx(k,:)) + (c_dbdy(j,:).*c_dbdy(k,:));
            % r = ((c_dbdx(:,j).*c_dbdx(:,k)) + (c_dbdy(:,j).*c_dbdy(:,k)))';
            
            % sum up the residuals for each node
            t_val(count) = r*c_weight;
            
            % count one up
            count = count + 1;
        end
    end
end
end
