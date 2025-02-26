function [full_sol, a_sel] = read_p2a(full_sol, opts)
%READ_A2P Convert the local sound speed to pressure assuming reference
%pressure and temperature

% specify the reference values
t_ref = opts.ref_temperature; 
p_ref = opts.ref_pressure; 
kappa = opts.kappa; 
R = opts.R;

% reshape the input field
n_vars = length(full_sol.variables);
n_nodes = size(reshape(full_sol.values(:,1), [], n_vars), 1);
a_sel = zeros(n_nodes, n_vars);
p_sel = ismember(full_sol.variables, {'p'});
a_sel(:, p_sel) = 1; a_sel = logical(reshape(a_sel, [], 1));

% assume slight compressibility of the fluid and convert pressure into
% temperature using a isentropic relation
is_rel = @(p) t_ref * ((p+p_ref)/p_ref).^((kappa-1)/kappa);
speed_of_sound = @(t) sqrt(kappa*R*t);

% decompose velocity and pressure field seperately
t_field = is_rel(full_sol.values(a_sel, :));
full_sol.values(a_sel, :) = speed_of_sound(t_field);

if any(p_sel)
    full_sol.variables{p_sel} = 'a'; 
else
    a_sel = zeros(n_nodes, n_vars);
    a_sel(:, ismember(full_sol.variables, {'a'})) = 1;
    a_sel = logical(reshape(a_sel, [], 1));
end
end

