function dir_sel = get_dir_points(fe_model)
%GET_DIR_POINTS Summary of this function goes here
%   Detailed explanation goes here
% get the start ids of the different variables
aid = find(ismember(fe_model.sys_variables, 'a')); 
ast = (aid - 1) * fe_model.mesh.n_node;
uid = find(ismember(fe_model.sys_variables, 'u')); 
ust = (uid - 1) * fe_model.mesh.n_node;
vid = find(ismember(fe_model.sys_variables, 'v')); 
vst = (vid - 1) * fe_model.mesh.n_node;

% select the nodes associated with the different dirichlet bounds
u_inlet_main = fe_model.mesh.lines(1).nodes + ust; % main inlet velocity
v_inlet_main = fe_model.mesh.lines(1).nodes + vst; % main inlet velocity
a_outlet = fe_model.mesh.lines(2).nodes + ast; % outlet pressure
v_inlet_hub = fe_model.mesh.lines(3).nodes + vst; % hub inlet velocity
v_inlet_shr = fe_model.mesh.lines(4).nodes + vst; % shr inlet velocity
dir_sel = [u_inlet_main, v_inlet_main, a_outlet, v_inlet_hub, v_inlet_shr];
end