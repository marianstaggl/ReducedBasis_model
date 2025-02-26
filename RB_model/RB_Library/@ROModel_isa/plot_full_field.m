function plot_full_field(self, prim_field, visc_field)
%PLOT_FULL_FIELD Plot the full field of the solution

% concatenate the field and the names
all_names = [self.fe_model.sys_variables(:)', self.visc_variables(:)'];
all_fields = [reshape(prim_field, [], 1); reshape(visc_field, [], 1)];
self.fe_model.mesh.plot_field(all_fields, all_names);
end

