function f_loc = get_fun_val(self, f_val, loc)
%GET_FUNCTION evaluate the function at the loc-points
%   the function specified by the vector f_val is evaluated at the loc
%   which is eighter the integration or the reference points

% get the function values at the reference points
f_ref = f_val(self.elems);

% if the location is the reference points
if isequal(loc, 'ref')
    f_loc = f_ref;
    
% if the location is on the integration points
elseif isequal(loc, 'int')
    % get the basefunction at the integration points
    int_points = self.ref_elem.int_points;
    b_f = self.ref_elem.get_b(int_points);
    
    % sum up all the base functions at the int points
    f_loc = f_ref * b_f;
else
    error('the specified query location is not available')
end
end

