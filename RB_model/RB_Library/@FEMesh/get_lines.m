function all_lines = get_lines( nodes, edges, elems )
%GET_LINES sort given edges into a line
%   the edges define the boundary of the mesh domain. they contain the node
%   ids of the nodes which form them. this function converts the edges
%   array into multiple lines which contain the corresponding nodes in
%   sorted manner.

% check if a line has 3 entries if it has two add a dummy column
if size(edges,2)==2, edges = [edges, zeros(size(edges,1),1)]; end

% find the available edge ids
ids = unique(edges(:,end)); 

% loop over all ids and find sorted lines
l_c = 1;
for i=1:length(ids)
    
    % define a template for the line
    c_l = struct('id', ids(i), 'nodes', [], 'pos', [],...
        'edges', [], 'elems', [], 'sub_id', [], 'n1', [], 'n2', []);
   
    % get a selector for the nodes
    sel = edges(:,end)==ids(i); edg_id = find(sel);
    
    % find the lines for the corresponding ids
    lines = get_single_line( c_l, nodes, edges(sel,1:end-1), edg_id );
    
    % save the lines
    all_lines(l_c:l_c + length(lines) - 1) = lines;
    
    % count up
    l_c = l_c + length(lines);
end

% associate an element to each line edge
all_lines = associate_elem(all_lines, edges(:,1:end-1), elems);

% add the wall normal to each line edge
all_lines = add_normals(all_lines, edges(:,1:end-1), nodes);
end

function all_l = add_normals(all_l, edges, node)
% loop over all the lines
for i=1:numel(all_l)
    % pick the current line to get the node ids
    c_edg = edges(all_l(i).edges, :);
    s_id = c_edg(:,1); f_id = c_edg(:,end);
    
    % get the n1 and n2 vecs
    all_l(i).n2 = normr(node(s_id,:) - node(f_id,:));
    all_l(i).n1 = normr([-all_l(i).n2(:,2), all_l(i).n2(:,1)]);
end
end

function all_lines = associate_elem(all_lines, edges, elems)
% associate each edge with an element

% loop over all the lines and their edges
for i=1:numel(all_lines)
    for j=1:length(all_lines(i).edges)
        % get the associated element
        edg_id = all_lines(i).edges(j);
        [row_1, ~] = find(elems == edges(edg_id,1));
        [row_2, ~] = find(elems == edges(edg_id,end));
        
        % find the corresponding edge
        elem_id = intersect(row_1, row_2);
        all_lines(i).elems = [all_lines(i).elems, elem_id];
        [~, sub] = ismember(edges(edg_id,:), elems(elem_id,:));
        all_lines(i).sub_id = [all_lines(i).sub_id, sub'];
        
        % check result
        if length(elem_id) ~= 1, error("this shouldn't be possible"); end
    end
end
end

function lines = get_single_line( lines, nodes, edges, edges_id )
% erzeugt zusammenh?ngende linien aus unsortierten kanten

% get the end nodes of each edge
end_nodes = edges(:,[1, end]); npe = size(edges,2);

% falls offene enden vorhanden sind wird dort angefangen
lines = repmat(lines,1,100); l_c = 1;
while ~isempty(end_nodes)
    
    % shortcuts anlegen
    id_new = reshape(end_nodes,1,[]); unique_new = unique(id_new);
    
    % nachsehen ob es freie enden gibt
    [n,~] = histc(id_new, unique_new); [~,n_idx] = sort(n);
    
    % startknoten festlegen
    curr_node = unique_new(n_idx(1));
    
    % so lange weitersuchen bis das ende der linie erreicht ist
    while true
        
        % abspeichern des endknoten
        lines(l_c).nodes = [lines(l_c).nodes, curr_node];
        
        % finden der momentanen zeile und des n?chsten knoten
        [new_row,new_col] = ind2sub(size(end_nodes),find(end_nodes == curr_node));
        
        % kontrolle ob man am ende angelangt ist
        if isempty(new_row), break; end
        
        % raussuchen des neuen knoten
        curr_node = end_nodes(new_row(1), ~(new_col(1)-1)+1);
        
        % abspeichern der inneren knoten
        ids = (new_col==1)*(2:1:npe-1) + (new_col==2)*(npe-1:-1:2);
        lines(l_c).nodes = [lines(l_c).nodes, edges(new_row(1), ids)];
        lines(l_c).edges = [lines(l_c).edges, edges_id(new_row(1))];
        
        % l?schen der zeile aus der conn_mat und counter raufz?hlen
        end_nodes(new_row(1),:) = []; edges(new_row(1),:) = []; 
        edges_id(new_row(1)) = [];
    end
    
    % set the nodes positions
    lines(l_c).pos = nodes(lines(l_c).nodes, :);
    
    % counter raufz?hlen
    l_c = l_c + 1;
end

% zurechtschneiden der linien
lines = lines(1:l_c-1);

end
