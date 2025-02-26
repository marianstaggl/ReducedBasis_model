function visc_state = init_viscos(self, target_field)
%INIT_VISCOS Initialize the viscos variables of the reduced system
%   the desired boundary conditions are approximated by scaling the
%   boundary modes of the ROModel. The scaling is done by least squares
%   fit.

% use a gappy pod approach to scale the boundary modes
bound_sel = true(size(self.visc_shape, 1), 1); %self.bound_nodes.viscos;

mean_mode = self.visc_shape(:, 1);
bound_mode = self.visc_shape(:, 2:sum(self.visc_bound));
res_field = target_field - mean_mode;

visc_state = bound_mode(bound_sel, :) \ res_field(bound_sel);
visc_state = [1; visc_state];
end

