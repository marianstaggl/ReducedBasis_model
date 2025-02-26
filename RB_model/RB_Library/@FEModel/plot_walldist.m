function plot_walldist(self)
%PLOT_FIELD plot the result of a calculation

% plot the given field on the current mesh
self.mesh.plot_field(self.wall_dist, {'wall distance'})
end

