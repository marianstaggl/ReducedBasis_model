function new_res = interp_res(par_tri, dataset, query_point)
%INTERPOLATE_RES Interpolate within the parameterspace

% enforce a row vector for the query point
query_point = reshape(query_point, 1, []);

% get the parameters together
par_mat = [dataset.geo_redpar]';

% get the polyeder enclosing the parameterpoint
[t,P] = tsearchn(par_mat, par_tri, query_point);

% check if there is a enclosing symplex
if isnan(t)
    % get the id of the closest point
    [c_pts,P] = closest(par_mat, par_tri, query_point);
else
    % get the ids of the enclosing simplex corners
    c_pts = par_tri(t(1),:);
end

% clone some values from the old datasets
new_res.name = 'interp_res';
new_res.geo_def = dataset(1).geo_def;
new_res.geo_base = dataset(1).geo_base;
new_res.geo_redpar = query_point';

% reconstruct the full parametrization
new_res.geo_parameter = new_res.geo_def + new_res.geo_base*new_res.geo_redpar;

% get all the fields of the dataset and the already implemented ones
all_field = fields(dataset); new_field = fields(new_res);

% get the differences in the set
c_field = setdiff(all_field, new_field);

% loop over all of the fields and interpolate the new result
for i=1:numel(c_field)
    c_base = [dataset.(c_field{i})];
    new_res.(c_field{i}) = c_base(:,c_pts)*P'; 
end

% order the field names of the new result
new_res = orderfields(new_res, dataset(1));
end

function [c_pts, P] = closest(P, T, PQ)

% search for the closest point
c_pts = dsearchn(P, T, PQ);

% set the weight to 1 at this point
P = 1;
end