function [ new_nodes, new_edges, new_elements ] = get_mesh( nodes, edges, elements, order )
%GET_MESH according to the choosen order the initial grid is split up
%   the standard element has only the four nodes at the corners but to use
%   higher order elements we need additional degrees of freedom. this is
%   achieved by inserting new nodes into the mesh

% create shortcuts for element number and the nodes per element
nel = size(elements,1); npelo = size(elements,2); np = size(nodes,1);
nple = (order + 1); npln = nple^2;

% sort the nodes of the elements matrix
elements = sort_element_corners(elements, nodes);

% initial check: every element should contain only 4 corners and one zone
% id making 5 entries per element in the elements matrix
if ~(npelo==5 || npelo==4)
    error('you can pass only linear elements with 4 corner nodes!');
end

% initial check if the order of the mesh should be 1
if order ==1
    new_nodes = nodes(:,1:end-1); new_edges = edges; new_elements = elements(:,1:end-1); return; 
end

% get the subids for the edges of the reference element
[edge_ids, corner_ids] = get_edges_and_corners( order );
edge_ids = [edge_ids; fliplr(edge_ids)];

% set the new_elements array
new_elements = nan(size(elements,1), npln);
new_elements(:,corner_ids) = elements(:,1:4);

% set the new_edges array
new_edges = nan(size(edges,1), sqrt(npln) + 1);
new_edges(:,1) = edges(:,1); new_edges(:,end-1:end) = edges(:,end-1:end);

% preallocate the new nodes array
new_nodes = [nodes(:,1:2); zeros(100000,2)];

% define a base element for the splitting of the 4 corner element
[n_opos, n_subid] = FElement.get_ref_points(order); 
temp_ref = FElement(1);

% loop over each element
for i=1:nel
    
    % get the current elements corner nodes
    c_e = elements(i,1:4); c_n = nodes(c_e,1:2);
    n_globid = nan(sqrt(npln), sqrt(npln));
    n_globsel = isnan(new_elements(i,:));
    
    % map the reference elements nodes onto the current one
    n_npos = temp_ref.map_pnts( c_n, n_opos );
    
    % get the new global node ids
    n_globid(n_globsel) = (1:1:sum(n_globsel)) + np;
    
    % add the new nodes to the element
    new_elements(i,n_globsel) = n_globid(n_globsel);
    
    % elements sharing an edge with the current one
    [n_elid, n_eedgid ,c_eedgid] = match_element_edges(i,...
        elements(i,1:4), elements(:,1:4));
    
    % loop over the neighboring elements
    for j=1:length(n_elid)
        % get the elements sub and global node ids
        c_sedg = n_subid(edge_ids(n_eedgid(j),:));
        c_gedg = n_globid(edge_ids(c_eedgid(j),:));
        
        % select only entries with nans
        sel_el = ~isnan(c_gedg);
        new_elements(n_elid(j),c_sedg(sel_el)) = c_gedg(sel_el);
    end
    
    
    % get the boundarys sharing an edge with the current element
    [n_bid, n_bedgid ,c_bedgid] = match_element_edges(nan, elements(i,1:4),...
        edges(:,1:2));
    
    % loop over the neighboring boundarys
    for j=1:length(n_bid)
        % get the elements sub and global node ids
        c_sedg = n_subid(edge_ids(n_bedgid(j),:));
        c_gbdg = n_globid(edge_ids(c_bedgid(j),:));
        
        % select only entries with nans
        sel_bo = ~isnan(c_gbdg);
        new_edges(n_bid(j), c_sedg(sel_bo)) = c_gbdg(sel_bo);
    end
    
    % at last add the nodes to the array
    new_nodes(n_globid(n_globsel),1:2) = n_npos(n_globsel,:);
    np = np + sum(n_globsel);
end

% cut the nodes
new_nodes = new_nodes(1:np,:);
end

function [n_elid, n_edgid ,c_edgid] = match_element_edges(ci, c_elem, elements)

% define the edges for the element with 4 nodes
c_edg = [1, 2; 2, 4; 3, 4; 1, 3; 2, 1; 4, 2; 4, 3; 3, 1];

% preallocate the arrays
n_elid = zeros(4,1); n_edgid = zeros(4,1); c_edgid = (1:1:4)';

% loop over all of the 4 edges
for i=1:4
    
    % get the indices of the first and second nodes
    [i1,j1] = ind2sub(size(elements), find(elements == c_elem(c_edg(i,1))));
    [i2,j2] = ind2sub(size(elements), find(elements == c_elem(c_edg(i,2))));
    
    % get the neighbor element id and the edge they share
    temp_id = intersect(i1,i2); temp_id =  temp_id(temp_id~=ci);
    
    % check if there are neighboring elements, otherwise continue
    if isempty(temp_id), continue; end
    n_elid(i) = temp_id;
    
    % get the nodes sub ids
    n_edg = [j1(i1 == n_elid(i)), j2(i2 == n_elid(i))];
    n_edgid(i) = find(ismember(c_edg,n_edg,'rows'));
end

% search for the zero elements
sel = ~(n_elid == 0);

% remove the 0 values
n_elid = n_elid(sel); n_edgid = n_edgid(sel); c_edgid = c_edgid(sel);
end

function elements = sort_element_corners(elements, nodes)

% sort the corner nodes of each 4point element
for i=1:size(elements,1)
    % get the current element and its corner nodes
    c_el = elements(i,1:4); c_no = nodes(c_el,1:2);
    
    % sort these nodes and save it into elements
    [~,sid] = FElement.get_sorted_corners(c_no);
    elements(i,1:4) = elements(i,sid);
end
end

function [edge_ids, corner_ids] = get_edges_and_corners( order )

% get the subids for the edges of the reference element
[~, subids] = FElement.get_ref_points( order );

% sort the subids into 4 edges
edge_ids = zeros(4,size(subids,1));
edge_ids(1,:) = subids(:,1);   edge_ids(2,:) = subids(end,:);
edge_ids(3,:) = subids(:,end); edge_ids(4,:) = subids(1,:);

% sort the subids into 4 corners
corner_ids = zeros(1,4);
corner_ids(1) = subids(1,1);   corner_ids(2) = subids(end,1);
corner_ids(3) = subids(1,end); corner_ids(4) = subids(end,end);
end

