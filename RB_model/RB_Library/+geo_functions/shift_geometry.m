function [geo_lines, shift] = shift_geometry(geo_lines)
%SHIFT_GEOMETRY

% get the mean zpos at the inlet
h_sel = ismember({geo_lines.name}, 'hub'); z_h = geo_lines(h_sel).pos(1,1);
s_sel = ismember({geo_lines.name}, 'shroud'); z_s = geo_lines(s_sel).pos(1,1);

% check if they are identical
if ~isequal(z_h, z_s), warning('z_h and z_s are different...'); end

% get the shift as mean value
shift = mean([z_h, z_s]);

% shift hub and shroud lines
geo_lines(h_sel).pos(:,1) = geo_lines(h_sel).pos(:,1) - shift;
geo_lines(s_sel).pos(:,1) = geo_lines(s_sel).pos(:,1) - shift;

% loop over the blade lines and shift them
b_id = find(contains({geo_lines.name}, 'blade'));
for i=1:length(b_id)
   geo_lines(b_id(i)).pos(:,3) =  geo_lines(b_id(i)).pos(:,3) - shift;
end
end

