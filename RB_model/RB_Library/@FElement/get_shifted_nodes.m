function nodes = get_shifted_nodes(old_supp, new_supp, nodes)
%GET_SHIFTED_NODES shift the nodes in a bilinear manner
%   the function takes the old corner nodes and the new corner nodes and
%   deformes the element accordingly. the other nodes are deformed
%   according to a bilinear transformation

% the shifting always works with bilinear shape functions
c_cmat = [ones(4,1), old_supp(:,1), old_supp(:,2), old_supp(:,1).*old_supp(:,2)];
c_smat = [ones(size(nodes,1),1), nodes(:,1), nodes(:,2), nodes(:,1).*nodes(:,2)];

% which coefficients are necessary to shift the corner nodes
shift_c = c_cmat\new_supp;

% shift the other nodes accordingly
nodes = c_smat*shift_c;

end

