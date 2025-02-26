function [new_res, int_flag] = interp_nD(par_tri, par_mat, data, query_points)
%INTERP_ND Interpolate within the parameterspace
%   use a triangulated surface to interpolate the results

% enforce a row vector for the query point
new_res = zeros(size(data,1), size(query_points,1));
int_flag = false(1, size(query_points,1));
for i=1:size(query_points,1)
    [new_res(:, i), int_flag(i)] = interpolate_single(par_tri, ...
        par_mat, data, query_points(i,:));
end
end

function [new_res, int_flag] = interpolate_single(par_tri, par_mat, data, query_point)
% get the enclosing simplex of the query point
[t, P] = tsearchn(par_mat, par_tri, query_point);

% check if there is an enclosing simplex
if isnan(t)
    % get the id of the closest case (nearest neighbor interpolation)
    int_flag = false;
    [c_pts, P] = closest(par_mat, par_tri, query_point);
else
    % get the ids of the enclosing simplex corners
    int_flag = true;
    c_pts = par_tri(t(1),:);
end

% interpolate the result
new_res = data(:,c_pts)*P';
end

function [c_pts, P] = closest(P, T, PQ)

% search for the closest point
c_pts = dsearchn(P, T, PQ);

% set the weight to 1 at this point
P = 1;
end