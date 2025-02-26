function plot_dataset(mesh, dataset, fields)
%PLOT_DATASET plot the given fields of the given result

% deform the mesh according to the parametrization
mesh = functions.get_geo('mesh', mesh, 'par', dataset.geo_parameter);

% get the given fields from the dataset
data = cell(1,numel(fields));
for i=1:numel(fields), data{i} = dataset.(fields{i}); end

% reshape the fields to plug them into the function
data = cell2mat(data);

% plot the results
mesh.plot_field(data, fields);

end

