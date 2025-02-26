function f_interp = get_fun_interp(self, field_val, query_pts)
%GET_FUN_INTERP interpolate the field on query points
%   use the mesh to interpolate the field on the query points

% get the triangulated mesh
tri = self.get_triangulation();
f_interp = interpolate_barycentric(tri, self.nodes, field_val, query_pts);
end

function interp_data = interpolate_barycentric(tri, data_pos, data_val, interp_pos)
%INTERP_ND Interpolate within the parameterspace
%   use a triangulated surface to interpolate the results

% run some checks if any of the arrays is a dl one
if isa(data_pos, 'dlarray'), data_pos = extractdata(data_pos); end
if isa(interp_pos, 'dlarray'), interp_pos = extractdata(interp_pos); end

% get the enclosing simplex
[t, P] = tsearchn(data_pos, tri, interp_pos);

% use the barycentric coordinates for interpolation
int_flag = ~isnan(t);
interp_data = zeros(size(interp_pos, 1), 1);

% convert the array into dl if necessary
if isa(data_val, 'dlarray'), interp_data = dlarray(interp_data); end
interp_data(int_flag) = sum(data_val(tri(t(int_flag), :)) .* P(int_flag, :), 2);

% use nearest neighbor for extrapolations
ext_flag = isnan(t);
% [~, min_id] = min(pdist2(data_pos, interp_pos(ext_flag, :)), [], 1);
interp_data(ext_flag) = nan;
end
