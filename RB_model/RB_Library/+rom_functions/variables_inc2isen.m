function [full_sol, p_sel] = variables_inc2isen(fe_model, full_sol)
%VARIABLES_INC2ISEN Summary of this function goes here
% get a selector for pressure field
p_sel = zeros(fe_model.mesh.n_node, length(fe_model.sys_variables));
p_sel(:, ismember(fe_model.sys_variables, {'a', 'p'})) = 1;
p_sel = logical(reshape(p_sel, [], 1));

% assume slight compressibility of the fluid and convert pressure into
% temperature using a isentropic relation
t_ref = 298; p_ref = 101300; kappa = 1.4; R = 287;
is_rel = @(p) t_ref * ((p+p_ref)/p_ref).^((kappa-1)/kappa);
speed_of_sound = @(t) sqrt(kappa*R*t);

% decompose velocity and pressure field seperately
t_field = is_rel(full_sol(p_sel, :));
full_sol(p_sel, :) = speed_of_sound(t_field);
end

