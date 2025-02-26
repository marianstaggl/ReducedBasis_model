% modify the inlet boundary conditions
function fe_model = mod_inlet_flow(fe_model, par_vec)
%MODIFY_MODEL modify the fem_model

% set the options for the calculation
inlet_line = fe_model.mesh.lines(1);
s_vec = FEMesh.get_s_vec(inlet_line); 
u_inf = 0.1; nue = 1.4041e-05;
[x_vel, y_vel] = rom_functions.calc_blasius_solution(par_vec, s_vec, u_inf, nue);

uid = find(ismember(fe_model.sys_variables, 'u')); 
ust = (uid - 1) * fe_model.mesh.n_node;
vid = find(ismember(fe_model.sys_variables, 'v')); 
vst = (vid - 1) * fe_model.mesh.n_node;

% update the boundary conditions for the current doe step
fe_model.dir_bounds(inlet_line.nodes + ust, 2) = x_vel;
fe_model.dir_bounds(inlet_line.nodes + vst, 2) = y_vel;
end

