function labels = get_labels(self)
%GET_SOL_LAB Returns the ids of the solution entries
%   the solution is always returned as a vector and this function returns
%   the ids of the nodes plus the associated variable in a vector with the
%   same sorting

% preallocate some variables
nn = self.mesh.n_node; 
ns = self.sys_variables;

% get the matrix with node ids
nid = num2cell(repmat((1:1:nn)', 1, numel(ns)));
var = repmat(ns, nn, 1);

% get the labels into a cell array
labels = cellfun(@(x1, x2) [x1,'_', num2str(x2)],...
    var, nid, 'UniformOutput', false);
labels = reshape(labels, [], 1);
end

