function data = red_parspace(data, th)
%RED_PARSPACE reduce the parameter space of the results

% get all the geometry vectors
geo_mat = [data.geo_parameter];

% define the first geo to be default (we define an origin)
geo_def = geo_mat(:,1); geo_mat = geo_mat(:,2:end) - geo_def;

% if a threshold is passed, use svd to remove redundancies
if ~isnan(th)
    % remove redundancies in the geometrys
    [geo_mat, S, ~] = svd(geo_mat,'econ'); s = S*ones(size(S,1),1);
    
    % remove redundante or unimportant geometric parameters
    sel = s > th; geo_mat = geo_mat(:,sel).*s(sel)';
end

% project the results onto thre new base
red_geo = pinv(geo_mat)*([data.geo_parameter] - geo_def);

% project the geometries onto the new parspace
for i=1:numel(data)
    
    % transform the old parametrization
    data(i).geo_redpar = red_geo(:,i);
    
    % save the current base of the parametrization
    data(i).geo_base = geo_mat; data(i).geo_def = geo_def;
end
end

