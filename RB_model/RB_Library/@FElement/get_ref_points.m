function [node_pos, node_ids, corner_ids] = get_ref_points( order )
%GET_REF_POINTS Return the position of the lagrangian support points
%   the shape functions of the elements are defined as lagrangian
%   polynomials. this function returns the support points of these
%   polynomials as well as the ids of the nodes and the corners

% create the coordinate vectors
eta_vec = linspace(-1, 1, order + 1);
zeta_vec = linspace(-1, 1, order + 1);

% mesh a grid of element points on reference element
[ZETA, ETA] = meshgrid(eta_vec, zeta_vec);

% use the same routine as in get_b to sort them
node_pos = [reshape(ETA,[],1), reshape(ZETA,[],1)];
node_ids = reshape(1:1:(order + 1)^2, order + 1, order + 1);
corner_ids = [node_ids(1,1), node_ids(end,1), node_ids(1,end), node_ids(end,end)];

end