function [nodes, edges, elements] = read_meshfile( self, mshfile )
%READ_MESHFILE this function runs the meshfile and returns nodes and elems
%  the meshfile is a matlabskript and is run in this function. then the
%  meshfile properties are sorted into nodes, edges and elements and
%  returned

% run the meshfile as it is a matlab skript
run(mshfile)

% set the nodes, edges and elements
nodes = msh.POS; nodes = [nodes(:,1:self.dim), nodes(:,end)];
edges = msh.LINES; edges = [edges(:,1:self.dim), edges(:,end)];
elements = msh.QUADS; elements = [elements(:,1:2^self.dim), elements(:,end)];

end