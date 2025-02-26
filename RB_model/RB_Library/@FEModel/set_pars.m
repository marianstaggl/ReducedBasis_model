function set_pars(self, geo_pars, bound_pars)
%SET_PARS Set a new parametrization of the model
%   a new parameter combination is passed and the mesh and the boundary
%   conditions of the model are adjusted accordingly

% reshape the parameter vectors
geo_pars = reshape(geo_pars,[],1); bound_pars = reshape(bound_pars,[],1);

% get the new nodes positons according to the parametrization
new_nodes = self.g_par_fun(self.g_zero, self.g_par_base, geo_pars);

% add the new nodes into the mesh and adjust the mapping
self.mesh.set_nodes(new_nodes); self.set_mapping();

% get the new boundary values according to the parametrization
new_bounds = [self.b_par_base(:,1), self.b_par_base(:,2:end)*bound_pars];

% add the new dirichlet bounds to the femodel
self.dir_bounds  = new_bounds; self.i_DoF = ~ logical(new_bounds(:,1));

% add the neumann boundary conditions
self.neu_bounds = zeros(size(new_bounds));

% set the new parametrization into the object
self.geo_pars = geo_pars;
self.bound_pars = bound_pars;
end

