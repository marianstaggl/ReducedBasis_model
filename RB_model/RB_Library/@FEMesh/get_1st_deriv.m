function [dfdx,dfdy] = get_1st_deriv(self, f_val, loc)
%GET_DERIVATIVES Get the derivatives of the funcition specified by f_val at 
% the integration or at the reference points (specified by loc)

% get the function values at the reference points
f_ref = f_val(self.elems);

% if the location is the reference points
if isequal(loc, 'ref')
    dfdx = sum_values(f_ref, self.dbdx_ref);
    dfdy = sum_values(f_ref, self.dbdy_ref);
    
% if the location is on the integration points
elseif isequal(loc, 'int')
    dfdx = sum_values(f_ref, self.dbdx_int);
    dfdy = sum_values(f_ref, self.dbdy_int);
else
    error('the specified query location is not available')
end
end

% perform the contraction of the matricies
function dpde = sum_values(f_ref, dbde)
dpde = zeros(size(f_ref,1), size(dbde,2));
for i=1:size(f_ref, 1)
    dpde(i,:) = dbde(:,:,i)'*f_ref(i,:)';
end
end
