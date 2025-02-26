function bound_nodes = get_bound_nodes(self, root_folder)
%GET_BOUND_NODES Read the file with the boundary conditions
%   a file is provided which specifies the edges, associated with a
%   boundary condition. the nodes are returned by the function.

% read the boundary file into a table
bound_table = read_boundaryfile(fullfile(root_folder, 'boundary.txt'));

% initialize the node selector and loop over the boundary entries
prim_bound = cell(1, size(bound_table, 1));
visc_bound = cell(1, size(bound_table, 1));
for i=1:size(bound_table,1)
    cvar = bound_table(i,:).var;
    clid = bound_table(i,:).id;
    psel = ismember(self.fe_model.sys_variables, cvar);
    vsel = ismember(self.visc_variables, cvar);

    % check if it is a primary or viscos variable
    sel = [self.fe_model.mesh.lines.id] == clid;
    all_nodes = {self.fe_model.mesh.lines(sel).nodes};
    all_nodes = cell2mat(all_nodes);

    if any(psel)
        cst = (find(psel) - 1) * self.fe_model.mesh.n_node;
        prim_bound{i} = all_nodes + cst;
    elseif any(vsel)
        cst = (find(vsel) - 1) * self.fe_model.mesh.n_node;
        visc_bound{i} = all_nodes + cst;
    end
end

% concatenate all the identified nodes
bound_nodes.primary = cell2mat(prim_bound);
bound_nodes.viscos = cell2mat(visc_bound);
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
    b_data = textscan(fid, '%f%s', 'CommentStyle', '%');
    fclose(fid);
catch me
    
    % when an error occurs close the file and then rethrow the error
    fclose(fid);
    rethrow(me);
end

% reshape the result of the boundary reading
b_data = horzcat(num2cell(b_data{1}), num2cell(b_data{2}));

% convert the dataset into a table
b_data = cell2table(b_data, 'VariableNames',{'id' 'var'});

end
