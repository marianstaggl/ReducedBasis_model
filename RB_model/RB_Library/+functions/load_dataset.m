function [data, mesh] = load_dataset(path, dataset, file_name, file_vars)
%LOAD_DATASET load the available datasets into a struct
% get the dataset folder
d_folder = fullfile(path, dataset);

% load the dataset into a struct array
data = load_results(d_folder, file_name, file_vars);

% add the parameter of the calculations
data = load_parameter(data, [d_folder '\geo_parameter.txt'], 12);

% create a mesh to plot the results
mesh = functions.get_geo();
end

function cfd_res = load_results(folder, file_name, file_vars)

% list all the folders of the dataset (containing sub)
folders = dir(fullfile(folder, 'sub*'));

% initialize the field values of the struct
f_vals = repmat({[]}, 1, numel(file_vars));

% initalize the struct
cfd_res = vertcat(file_vars, f_vals); cfd_res = struct(cfd_res{:}); cfd_res.name = '';

% preallocate it with the necessary number
cfd_res = repmat(cfd_res, 1, numel(folders));

% loop over all of the folders to get the results
for i=1:numel(folders)
    
    % get the resultfile name
    filename = fullfile(folders(i).folder, folders(i).name, file_name);
    
    % open the file in a try catch block
    try
        % open the file and read the boundary data
        fid = fopen(filename,'r');
        data = textscan(fid,repmat('%f', 1, numel(file_vars)), 2145,...
            'HeaderLines',2,'Delimiter','\t','CollectOutput',1,'EndOfLine','\r\n');
        fclose(fid);
    catch me
        
        % when an error occurs close the file and then rethrow the error
        fclose(fid);
        rethrow(me);
    end
    
    % convert it into a matrix
    data = data{1};
    
    % put the results into a struct
    for j=1:numel(file_vars), cfd_res(i).(file_vars{j}) = data(:,j); end
    
    % and set the name
    cfd_res(i).name = folders(i).name;  
%     cfd_res(i).x = data(:,1); cfd_res(i).y = data(:,2); cfd_res(i).z = data(:,3);
%     cfd_res(i).rho = data(:,4); cfd_res(i).u = data(:,5); cfd_res(i).v = data(:,6);
%     cfd_res(i).w = data(:,7); cfd_res(i).p = data(:,8); 
end
end

function data = load_datablock(fid)
% read the datablock of a given file. the datablock is assumed to be
% numeric, one continuous part and of table-like structure. the start and 
% endline is identified automatically

% preallocate a cellarray to save the file data and the counter
data = cell(1,1e6); i = 1;

% loop over the lines of the file
while ~feof(fid), data{i} = strsplit(fgetl(fid)); i=i+1; end

% get the length of the lines
c_len = cellfun(@length, data); m_len = median(c_len(c_len~=0));

% find the biggest table-like block
ids = reshape(find(diff([0;(c_len == m_len)';0])~=0),2,[]);

% pick the identified datablock
data_block = {data{ids(1):ids(2)-1}}; data_block = vertcat(data_block{:});

% convert the entries to numeric type
data = cell2mat(cellfun(@ str2num, data_block, 'UniformOutput', false));
end

function cfd_res = load_parameter(cfd_res, filename, no_par)

% open the file in a try catch block
try
    % open the file and read the boundary data
    fid = fopen(filename,'r');
    data = textscan(fid, ['%s', repmat('%f', 1, no_par)],...
        'EndOfLine','\r\n','CollectOutput',1);
    fclose(fid);
catch me
    
    % when an error occurs close the file and then rethrow the error
    fclose(fid);
    rethrow(me);
end

% loop over all the cfd results and add the parametrization
for i=1:numel(cfd_res)
    
    % find the corresponding line in the parameterfile
    sel = ismember(data{1}, cfd_res(i).name);
    
    % check if there is exactly one match
    if sum(sel) ~= 1, error("parameterfile doesn't match..."); end
    
    % add the geometry parameters to the datafile
    cfd_res(i).geo_parameter = reshape(data{2}(sel,:),[],1);
end
end