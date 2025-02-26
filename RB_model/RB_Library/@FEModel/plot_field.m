function plot_field(self, field, varargin)
%PLOT_FIELD plot the result of a calculation

% plot the given field on the current mesh
self.mesh.plot_field(field, self.sys_variables, varargin{:})
end

