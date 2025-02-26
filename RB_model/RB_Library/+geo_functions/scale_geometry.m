function  [geo_lines, scale] = scale_geometry(geo_lines)
%NORMALIZE_GEO scale the geometry according to mean inlet radius

% get the inlet radius of hub and shroud
r_h = geo_lines(ismember({geo_lines.name}, 'hub')).pos(2,2);
r_s = geo_lines(ismember({geo_lines.name}, 'shroud')).pos(2,2);

% get the mean radius and get the scale
r_m = mean([r_h, r_s]);
scale = 0.25/r_m; % the mean radius is converted to 0.25

% scale all the nodes positions of the lines
for i=1:numel(geo_lines), geo_lines(i).pos = scale*geo_lines(i).pos; end
end

