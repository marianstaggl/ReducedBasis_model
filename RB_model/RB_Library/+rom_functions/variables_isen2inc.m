function [full_sol, p_sel] = variables_isen2inc(fe_model, full_sol)
%VARIABLES_ISEN2INC Summary of this function goes here
% get a selector for pressure field
p_sel = zeros(fe_model.mesh.n_node, length(fe_model.sys_variables));
p_sel(:, ismember(fe_model.sys_variables, {'a', 'p'})) = 1;
p_sel = logical(reshape(p_sel, [], 1));

% assume slight compressibility of the fluid and convert pressure into
% temperature using a isentropic relation
t_ref = 298; p_ref = 101300; kappa = 1.4; R = 287;
temp_rel = @(a) a.^2 / (kappa*R) ;
is_rel = @(t) p_ref * (t/t_ref).^(kappa/(kappa-1)) - p_ref;

% decompose velocity and pressure field seperately
t_field = temp_rel(full_sol(p_sel, :));
full_sol(p_sel, :) = is_rel(t_field);
end

