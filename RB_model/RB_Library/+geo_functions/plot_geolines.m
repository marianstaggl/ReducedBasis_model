function h = plot_geolines(geo_lines, varargin)
%PLOT_GEOLINES plot the geometry lines

% make a figure and format
h = figure(); ax1 = subplot(2,1,1); hold on; grid on; axis equal;
ax2 = subplot(2,1,2); hold on; grid on; axis equal;

% get the hub line and plot it
hub_line = geo_lines(ismember({geo_lines.name}, 'hub'));
shr_line = geo_lines(ismember({geo_lines.name}, 'shroud'));

% plot hub and shroud lines
subplot(2,1,1); plot(hub_line.pos(:,1), hub_line.pos(:,2), '.');
subplot(2,1,1); plot(shr_line.pos(:,1), shr_line.pos(:,2), '.');

% plot the contour lines of the geo
l_n = {geo_lines.name}; b_lines = geo_lines(contains(l_n, 'blade'));
for i=1:numel(b_lines)
   subplot(2,1,1); plot(b_lines(i).pos(:,3), b_lines(i).pos(:,1))
   subplot(2,1,2); plot(b_lines(i).pos(:,3), b_lines(i).pos(:,2), '.');
end

% link the axes of the subplots
linkaxes([ax1, ax2], 'x');

% if needed save the plot as figure
if ~isempty(varargin), savefig(h, varargin{1}); end
end

