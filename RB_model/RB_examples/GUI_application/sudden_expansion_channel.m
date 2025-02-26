%% clear the workspace and define parameters
% clear the workspace and add the library
clc; clear all; close all;
addpath('..\..\RB_Library')
root_path = '..\..\RB_data\rom_root_sudden_expansion_channel';

%% setup the reduced model and load the data
% initialize the reduced order model
ro_model = ROModel_isa(root_path);
ro_model.set_base_fun('em19_gm4');
ro_model.set_sys_mat();
ro_model.set_deim_elem();

% set up the geo function
geo_base = load('+util_fun\geo_base.mat').geo_base;
geo_fun = @(x) reshape(geo_base(:,1) + x*geo_base(:,2), [], 2);

%% start the GUI
par_list = {'geo_par', 'vel_par', 'visc_par'};
eval_fun = @(ro_model, par_struct) run_calculation(ro_model,...
    par_struct, geo_fun);
ROMApplication(ro_model, eval_fun, par_list); 

%% functions block
function s = run_calculation(ro_model, par_struct, geo_fun)
geo_par = par_struct.geo_par;
vel_par = par_struct.vel_par;
visc_par = par_struct.visc_par;

% scale the parameters accordingly
geo_limits = [-25, 5]; vel_limits = [5e-4, 1e-2]; visc_limits = [1, 20];
geo_par = min(geo_limits) + geo_par * abs(diff(geo_limits));
vel_par = min(vel_limits) + vel_par * abs(diff(vel_limits));
visc_par = min(visc_limits) + visc_par * abs(diff(visc_limits));

% deform the geometry of the fe model
ro_model.fe_model.mesh.set_nodes(geo_fun(geo_par));
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