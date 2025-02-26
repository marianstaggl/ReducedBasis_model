function [d2fdx2,d2fdy2] = get_2nd_deriv(self,f_val, loc)
%GET_DERIVATIVES Get the derivatives of the funcition specified by f_val at 
% the integration or at the reference points (specified by loc)

% get the function values at the reference points
f_ref = f_val(self.elems);
    
% if the location is on the integration points
if isequal(loc, 'int')
    d2fdx2 = sum_values(f_ref, self.d2bdx2_int);
    d2fdy2 = sum_values(f_ref, self.d2bdy2_int);
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

