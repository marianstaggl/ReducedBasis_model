function plot_reduced_field(self, prim_red, visc_red)
%PLOT_REDUCED_FIELD Plot the reduced solution

% project the reduced solutions onto the modes
prim_field = self.prim_shape * reshape(prim_red, [], 1);
visc_field = self.visc_shape * reshape(visc_red, [], 1);

% plot the full field
self.plot_full_field(prim_field, visc_field);
end

