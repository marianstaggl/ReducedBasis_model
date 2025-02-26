function int_val = get_integral(self, f_val)
%GET_INTEGRAL integrate the function defined by f_val

% integrate the function with the weights
fval_int = self.get_fun_val(f_val, 'int');
int_val = sum(fval_int.*self.weights, 'all');
end

