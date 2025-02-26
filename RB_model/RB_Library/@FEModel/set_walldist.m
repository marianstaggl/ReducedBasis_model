function set_walldist(self, wall_ids)
%SET_WALLDIST set the wall distance property of the mesh
%   for each node, the distance to the closest wall is calculated

% collect all the nodes associated to a wall
line_sel = ismember([self.mesh.lines.id], wall_ids);
wall_ids = [self.mesh.lines(line_sel).nodes];

% calculate the distance of each node to the wall nodes
wall_nodes = self.mesh.nodes(wall_ids, 1:2);
all_nodes = self.mesh.nodes(:, 1:2);
D = pdist2(all_nodes, wall_nodes, 'euclidean');

% save the wall distance within the object
[self.wall_dist, self.wall_ids] = min(D,[],2);
end

function [dbdx, dbdy] = map_elems(edg_elem, edg, node, ref_elem, par_elem)
    % sub function to map the derivatives of one single element
    function  [dbdx_edg, dbdy_edg] = map_single( dpdx_ref, dpdy_ref, dbde_ref, dbdz_ref, support_nodes )
        % map the properties of a single element
        
        % preallocate the arrays
        dbdx_edg = zeros(size(dpdx_ref)); 
        dbdy_edg = zeros(size(dpdx_ref));
        
        % loop over all of the integration points
        for j=1:size(dpdx_ref, 2)
            % get the inverse jacobian matrix in reference points
            ijac_ref = inv((support_nodes'*[dpdx_ref(:,j), dpdy_ref(:,j)])');
            
            % map the derivatives of the reference points
            dbdx_edg(:,j) = ijac_ref(1,1)*dbde_ref(:,j) + ijac_ref(1,2)*dbdz_ref(:,j);
            dbdy_edg(:,j) = ijac_ref(2,1)*dbde_ref(:,j) + ijac_ref(2,2)*dbdz_ref(:,j);
        end
    end

% preallocate some variables
n_base = size(edg_elem, 2);
n_pts = size(edg, 2);
n_elem = size(edg_elem, 1);

% preallocate some of the arrays of dim (n_base x n_points x n_elem)
dbdx = zeros(n_base, n_pts, n_elem);
dbdy = zeros(n_base, n_pts, n_elem);

% get the basefunctions on the edg points
dbde_edg = ref_elem.get_dbde(ref_elem.ref_points, 1);
dbdz_edg = ref_elem.get_dbdz(ref_elem.ref_points, 1);
dpdx_edg = par_elem.get_dbde(par_elem.ref_points, 1);
dpdy_edg = par_elem.get_dbdz(par_elem.ref_points, 1);

% loop over the edges and elements
for i=1:size(edg,1)
    % get the subids of the edg nodes within the element
    [~, sub] = ismember(edg(i,:), edg_elem(i,:));
    
    % get the support nodes of the element
    support = node(edg_elem(1,:), :);
    [dbdx(:,:,i), dbdy(:,:,i)] = map_single(dpdx_edg(:,sub),...
        dpdy_edg(:,sub), dbde_edg(:,sub), dbdz_edg(:,sub), support);
end
end

function [n1_vecs, n2_vecs, lens] = get_normals(edg, node)
% pick the current line to get the node ids
s_id = edg(:,1); f_id = edg(:,end);
lens = sqrt(sum((node(s_id,:) - node(f_id,:)).^2, 2));

% calculate the line normals (at the center of the pieces)
n2_vecs = normr(node(s_id,:) - node(f_id,:)); % p2p vecs
n1_vecs = normr([-n2_vecs(:,2), n2_vecs(:,1)]); % normalized normals
end

function edg_elem = associate_elem(elem, edg)

% preallocate the edge element array
edg_elem = zeros(size(edg,1), size(elem,2));

% loop over the wall id nodes and save the elements
for i=1:size(edg,1)
    % get the associated element
    [row_1, ~] = find(elem == edg(i,1));
    [row_2, ~] = find(elem == edg(i,end));
    
    % find the corresponding edge
    elem_id = intersect(row_1, row_2);
    edg_elem(i,:) = elem(elem_id, :);
end
end

function edg = split_line(line, order)
% split the line into edges
edg = reshape(line.nodes(1:end-1), order, [])';
edg = [edg, [edg(2:end,1); line.nodes(end)]];
end

