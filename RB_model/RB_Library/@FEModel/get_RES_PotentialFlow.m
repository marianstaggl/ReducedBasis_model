function Res = get_RES_PotentialFlow(self, curr_cond, ~)
%GET_OP_HEATTRANSF get the operator for the heat transfer equation

% check the systems variables (this only works for t)
if ~isequal(self.sys_variables, {'phi'})
    error('this function only works for a temperature field!')
end

% get the basefunction at the integration points
elems = self.mesh.elems;
dbdx = self.dbdx; dbdy = self.dbdy; 
weights = self.weights;

% get the current conditions and reshape them into a new order
phivec = reshape(curr_cond, [], 1);
phiref = phivec(self.mesh.elems); % get the potential in the reference points
dphidx = sum(permute(dbdx, [3, 2, 1]) .* permute(phiref, [1, 3, 2]), 3);
dphidy = sum(permute(dbdy, [3, 2, 1]) .* permute(phiref, [1, 3, 2]), 3);

% get the residuals in the integration points
phi_res = dphidx.*permute(dbdx, [3, 2, 1]) + ...
    dphidy.*permute(dbdy, [3, 2, 1]);
res_int = reshape(sum(phi_res.*weights, 2), 1, []);

% get the indices for the residual vector
k_idx = int64(reshape(elems, [], 1))';
Res = accumarray(k_idx', res_int');
end

% %%
% 
% % preallocate the matrix and set up the numbers
% [n_elem, n_base] = size(elems);
% 
% % preallocate the vectors for the sparse matrix assembly
% t_vec = curr_cond;
% temp = zeros(1,n_elem*n_base);
% t_val = temp; k_idx = temp;
% 
% % loop over all the elements
% count = 1;
% for i=1:n_elem
% 
%     % get the current element and its mappings
%     curr_elem = elems(i,:); 
% 
%     % get shortcuts to the mapped derivatives
%     c_weight = weights(i,:)'; c_dbdx = dbdx(:,:,i); c_dbdy = dbdy(:,:,i); 
% 
%     % loop over all of the elements nodes twice
%     for k=1:n_base
% 
%         r = (dphidx(i,:).*c_dbdx(k,:)) + (dphidy(i,:).*c_dbdy(k,:));
%         t_val(count) = r*c_weight;
%         k_idx(count) = curr_elem(k);
%         count = count + 1;
%     end
% end
% 
% %test_stuff
% Lres = accumarray(int64(k_idx'), t_val');
