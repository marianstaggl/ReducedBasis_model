function [var_values, var_names, all_param, coord_values] = ...
    read_ASCII(self, folder, var_assignment, varargin)
%READ_ASCII Read all ansys files and map them onto the mesh
%   the parametrization of the given results is specified by the filename.
%   the format is shown below. This file stems from a dataset which has two
%   parameters named ParName1, ParName2 with the values ParName1 = 0.01 and
%   ParName2 = 1.3.
%
%   #ParName1=0_01#ParName2=1_3.txt

ip = inputParser();
ip.addParameter('node_mapping', 'none');
ip.parse(varargin{:});

% read the nodemap if it is provided
node_map = get_node_map(self.fe_model, ip.Results.node_mapping);
[all_files, all_param] = get_sorted_files(folder);

% read the first file to get the variable names
[~, var_names, ~] = read_current_file(fullfile(folder,...
    all_files(1).name), var_assignment);

% read the ansys results
coord_values = cell(1, numel(all_files));
var_values = cell(1, numel(all_files));
for i=1:numel(all_files)
    % read the current file
    c_file = fullfile(folder, all_files(i).name);
    [coords, new_names, vals] = read_current_file(c_file, var_assignment);
    [~, idx] = ismember(var_names, new_names);
    vals = vals(:, idx);
    
    % get the values at the fe_mesh node positions
    if isempty(node_map)
        all_vars = functions.interp_invD(coords, vals,...
            self.fe_model.mesh.nodes, 1, 1);
    else
        all_vars = vals(node_map(:,2), :);
        coords = coords(node_map(:,2), :);
    end

% return the interpolated flow field
coord_values{i} = coords;
var_values{i} = all_vars;
end
end

% get the nodemap for the current file
function node_map = get_node_map(fe_model, node_map_file)
fe_nodes = 1:1:fe_model.mesh.n_node;

if ~isequal(node_map_file, 'none')
    node_map = sortrows(readmatrix(node_map_file));
    
    % ceck the consistency of the nodemap first
    if ~isequal(node_map(:,1)', fe_nodes)
        error(['there is a problem with the node map. every mesh '...
            'node needs an assigment to a node within the file.'])
    end
else
    node_map = [];
end
end

% read the given file and return the values
function [coords, other_name, vals] = read_current_file(filename,...
    var_assignment)
% specify variable names
file_vars = read_file_vars(filename, var_assignment);
res = DataObj.reader(filename, 'singletable', 'column_names', file_vars);

% get the names of the columns which should be extracted
coord_name = {'x', 'y'};
other_name = setdiff(fields(var_assignment), {'x', 'y'});

% loop over the system variables
coords = res.get_column(coord_name);
vals = res.get_column(other_name);
end

% read the variable names of the FLUENT-Files
function file_vars = read_file_vars(file_name, var_assignment)
% read the first line and get the variable names
try
    fid = fopen(file_name);
    file_vars = split(fgetl(fid));
    fclose(fid);
catch, fclose(fid); 
end

% remove empty cells from the field names
file_vars = reshape(file_vars(~cellfun('isempty', file_vars)), 1, []);

% loop over the fields in the assignment struct and replace the names
field_names = fieldnames(var_assignment);
for i=1:numel(field_names)
    selector = ismember(file_vars, var_assignment.(field_names{i}));
    file_vars(selector) = field_names(i);
end
end

% get all the filenames and parse the parameter values
function [all_files, all_param] = get_sorted_files(root_folder)
% get all the text files within the root folder
all_files = dir(fullfile(root_folder, '*.txt'));

% split up the names to get the parameters
all_param = [];
for i=1:length(all_files)
    all_param = [all_param; parse_filename(all_files(i).name)];
end
end

% translate the filename into the parametrization of the file
function c_param = parse_filename(filename)
% preallocate the table and split the filename accordingly
c_param = table();
c_name = strsplit(filename, {'#', '.'});
c_name = cellfun(@(x) strsplit(x, '='), c_name(2:end-1),...
    'UniformOutput', false);

% loop over the parameters and store them
for j=1:length(c_name)
    par_name = c_name{j}{1};
    par_val = str2double(strrep(c_name{j}{2}, '_', '.'));
    c_param.(par_name) = par_val;
end
end