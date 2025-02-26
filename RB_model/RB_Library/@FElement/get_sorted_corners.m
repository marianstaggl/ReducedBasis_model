function [nodes, sort_idx] = get_sorted_corners(nodes)
%GET_SORTED_CORNERS sorts the corner nodes in according to the element base
%   the corners need to be sorted to perform the mapping. the sequence of
%   the corners needs to match the seqence of the base functions. this
%   function removes the mean position of the element and sorts all nodes
%   according to the sketch below.
%
%   3_ _ _ _ _ _4
%   |           |
%   |           |
%   |           |
%   |           |
%   |_ _ _ _ _ _|
%   1           2

% set the center to 0
pos = nodes - mean(nodes,1);

% perform a coordinate transformation to polar coordinates
[theta,~] = cart2pol(pos(:,1),pos(:,2));

% sort the angles of the new system
[~, sort_idx] = sort(mod(theta,2*pi));

% sort the corners in the right sequence
sort_idx = sort_idx([3,4,2,1]);
nodes = nodes(sort_idx,:);

end

