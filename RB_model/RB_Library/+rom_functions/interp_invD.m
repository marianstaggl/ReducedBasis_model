function new_res = interp_invD(par_mat,data_mat, query_points, n_pts, w_exp)
%INTERP_INVD use inverse distance weighting to interpolate

% first calculate the distance of the sampling points to the support points
dist_mat = pdist2(par_mat, query_points);
new_res = zeros(size(data_mat,1), size(query_points,1));

% loop over the query points
for i=1:size(query_points, 1)
    % get the closest cases for each query point
    [val, idx] = sort(dist_mat(:, i));

    % select the first n_pts and calculate the weights
    inv_dist = (1 ./ val(1:n_pts)).^w_exp;
    weights = inv_dist / sum(inv_dist);
    c_idx = idx(1:n_pts);

    % calculate the new solution as a weighted sum
    new_res(:, i) = data_mat(:, c_idx) * weights;
end
end

