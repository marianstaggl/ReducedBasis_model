function s = evaluate_model(ro_model, par_vec, geo_deform)
% clip the parameters within 0 and 1
vel_par = functions.dl_clip(par_vec(1), 0, 1);
geo_par = functions.dl_clip(par_vec(2), 0, 1);
visc_par = functions.dl_clip(par_vec(3), 0, 1);

% scale the parameters accordingly
geo_limits = [-25, 5]; vel_limits = [5e-4, 1e-2]; visc_limits = [1, 20];
geo_par = min(geo_limits) + geo_par * abs(diff(geo_limits));
vel_par = min(vel_limits) + vel_par * abs(diff(vel_limits));
visc_par = min(visc_limits) + visc_par * abs(diff(visc_limits));

% deform the geometry of the fe model
ro_model.fe_model.mesh.set_nodes(geo_deform(geo_par));
ro_model.fe_model.set_mapping_red(ro_model.get_deform_red());

% initialize the outlet conditions
primary_field = dlarray(zeros(size(ro_model.prim_shape, 1), 1));
a_sel = find(ro_model.fe_model.get_var_sel({'a'}));
a_bound = a_sel(ro_model.fe_model.mesh.lines(4).nodes);
primary_field(a_bound) = 0.0080;

% get the inlet conditions for the rom model
u_sel = find(ro_model.fe_model.get_var_sel({'u'}));
u_bound = u_sel(ro_model.fe_model.mesh.lines(3).nodes);
primary_field(u_bound) = vel_par;
init_primary = ro_model.init_primary(primary_field);

% run the simulation of the reduced system
ro_opts = ro_model.get_opts('plot', false, 'shape2test', 3);
red_state = ro_model.run_simulation_dl(init_primary,...
    visc_par, ro_opts, 20);
s = ro_model.prim_shape * red_state;
end