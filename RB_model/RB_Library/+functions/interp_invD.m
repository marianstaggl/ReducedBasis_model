function new_res = interp_invD(data_pts, data_val, query_pts, n_pts, w_exp)
%INTERP_INVD use inverse distance weighting to interpolate

% first calculate the distance of the sampling points to the support points
dist_mat = pdist2(data_pts, query_pts);
new_res = zeros(size(query_pts,1), size(data_val,2));

% loop over the query points
for i=1:size(query_pts, 1)
    % get the closest cases for each query point
    [val, idx] = mink(dist_mat(:, i), n_pts);

    % select the first n_pts and calculate the weights
    inv_dist = (1 ./ val).^w_exp;
    weights = inv_dist / sum(inv_dist);

    % calculate the new solution as a weighted sum
    new_res(i, :) = data_val(idx, :)' * weights;
end
end