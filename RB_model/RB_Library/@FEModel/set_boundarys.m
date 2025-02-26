function set_boundarys( self, varargin )
%SET_BOUNDARYS set function to update the boundary data
%   the boundary data is updated via this function. you pass a boundary
%   file and the dir_bounds and the neu_bounds properties of the femodel
%   are updated accordingly. the boundary file has to have the following
%   format (to comment out a line use %):
%   edgeid(number)    type(n/d)    variable(u,v,p)    value

% check if the boundaryfile input is empty
if ~isempty(varargin)
    bound_data = read_boundaryfile(varargin{1});
elseif ~isempty(self.boundary)
    bound_data = self.boundary;
else, error('no boundary data specified'); end
bound_data_len = size(bound_data,1);

% get the number of nodes, variables and boundary ids
all_var = unique(bound_data.var);
num_var = numel(all_var);

% set the systems variables
self.sys_variables = all_var';

% preallocate the arrays
dir_bounds = zeros(self.mesh.n_node*num_var, 2); 
neu_bounds = zeros(self.mesh.n_node*num_var, 2);
rob_bounds = zeros(self.mesh.n_node*num_var, 2);

% loop over all of the boundarys
for i=1:bound_data_len
    
    % get the corresponding nodes
    sel = [self.mesh.lines.id] == bound_data.id(i);
    nodes = [self.mesh.lines(sel).nodes];
    
    % account for a number of variables > 1
    nodes = nodes + (find(ismember(all_var, bound_data.var(i)))-1)*self.mesh.n_node;
    
    % save it if is a neuman, dirichlet or robin bound
    dir_bounds(nodes,1) = dir_bounds(nodes,1) + isequal(bound_data.type(i),{'d'});
    neu_bounds(nodes,1) = neu_bounds(nodes,1) + isequal(bound_data.type(i),{'n'});
    rob_bounds(nodes,1) = rob_bounds(nodes,1) + isequal(bound_data.type(i),{'r'});
    
    % save the boundary values
    dir_bounds(nodes,2) = dir_bounds(nodes,2) + bound_data.value(i) * isequal(bound_data.type(i),{'d'});
    neu_bounds(nodes,2) = neu_bounds(nodes,2) + bound_data.value(i) * isequal(bound_data.type(i),{'n'});
    rob_bounds(nodes,2) = rob_bounds(nodes,2) + bound_data.value(i) * isequal(bound_data.type(i),{'r'});
end

% check if the dirichlet boundary conditions are intersecting
if any(dir_bounds(:,1)>1), warning('intersecting dirichlet boundary conditions (probably at the corners)'), end
if any(neu_bounds(:,1)>1), warning('intersecting neumann boundary conditions (probably at the corners)'), end
if any(rob_bounds(:,1)>1), warning('intersecting robin boundary conditions (probably at the corners)'), end

% write the vectors into the model
self.boundary = bound_data;
self.dir_bounds = dir_bounds; self.dir_bounds(:,1) = self.dir_bounds(:,1)>0;
self.neu_bounds = neu_bounds; self.neu_bounds(:,1) = self.neu_bounds(:,1)>0;
self.rob_bounds = rob_bounds; self.rob_bounds(:,1) = self.rob_bounds(:,1)>0;
self.i_DoF = ~ logical(dir_bounds(:,1));

end

function b_data = read_boundaryfile( boundaryfile )
%READ_BOUNDARYDATA read the given boundary file into a cellarray
%   the whole reading of the text file is embedded in a try catch block to
%   properly close the file if an error occurs. otherwise there are always
%   tons of files open without any handle to them

% use a try catch block to avoid opened files after error
try
    % open the file and read the boundary data
    fid = fopen(boundaryfile,'r');
    b_data = textscan(fid, '%f%c%s%f', 'CommentStyle', '%');
    fclose(fid);
catch me
    
    % when an error occurs close the file and then rethrow the error
    fclose(fid);
    rethrow(me);
end

% reshape the result of the boundary reading
b_data = horzcat(num2cell(b_data{1}), num2cell(b_data{2}),...
    num2cell(b_data{3}), num2cell(b_data{4}));

% convert the dataset into a table
b_data = cell2table(b_data, 'VariableNames',{'id' 'type' 'var' 'value'});

end