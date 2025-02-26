function [s_base, s_DoF, t_base, t_DoF] = get_rom_base(self, opts)
%GET_ROM_BASE Summary of this function goes here
%   Detailed explanation goes here

% cut the system matrices to the right size
n_shape = size(self.prim_shape, 2);
n_test = min([floor(opts.shape2test * n_shape), size(self.prim_test, 2)]);

% get the restricted base functions
s_base = self.prim_shape;
t_base = self.prim_test(:,1:n_test);

% specify the degrees of freedom
s_DoF = true(1, n_shape); s_DoF(1:sum(self.prim_bound)) = false;
t_DoF = true(1, n_test); t_DoF(1:sum(self.prim_bound)) = false;
end

