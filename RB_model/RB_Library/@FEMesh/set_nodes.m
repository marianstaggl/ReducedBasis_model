function set_nodes(self, new_nodes)
%SET_NODES new nodes can be inserted into the mesh

% first check if the number of nodes remains constant
if ~isequal(size(self.nodes), size(new_nodes))
    error('you can only change the positions of the nodes not their number!')
end

% then set the new nodes positions
self.nodes = new_nodes;

for i=1:numel(self.lines)
    self.lines(i).pos = self.nodes(self.lines(i).nodes, :);
end
end

