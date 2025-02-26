function prim_state = init_primary(self, target_field)
%INIT_PRIMARY Initialize the primary variables of the reduced system
%   the desired boundary conditions are approximated by scaling the
%   boundary modes of the ROModel. The scaling is done by least squares
%   fit.

% use a gappy pod approach to scale the boundary modes
bound_sel = self.bound_nodes.primary;

mean_mode = self.prim_shape(:, 1);
bound_mode = self.prim_shape(:, 2:sum(self.prim_bound));
res_field = target_field - mean_mode;

bound_init = pinv(bound_mode(bound_sel, :)) * res_field(bound_sel); % bound_mode(bound_sel, :) \ res_field(bound_sel);
domain_init = zeros(self.n_shape, 1); % fill the remainig initializations with zeros

% build the full initial state of the system
prim_state = [1; bound_init; domain_init];
end

