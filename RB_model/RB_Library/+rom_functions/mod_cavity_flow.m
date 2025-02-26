function fe_model = mod_cavity_flow(fe_model, par_vec)
%MODIFY_MODEL modify the fem_model

% change the inlet velocity through the cavity
cavity_line = fe_model.mesh.lines(5);
vid = find(ismember(fe_model.sys_variables, 'v'));
vst = (vid - 1) * fe_model.mesh.n_node;
fe_model.dir_bounds(cavity_line.nodes + vst, 2) = par_vec;
end