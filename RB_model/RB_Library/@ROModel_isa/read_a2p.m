function [full_sol, p_sel] = read_a2p(full_sol, opts)
%READ_A2P Convert the pressure to local sound speed assuming reference
%pressure and temperature

% specify the reference values
t_ref = opts.ref_temperature; 
p_ref = opts.ref_pressure; 
kappa = opts.kappa; 
R = opts.R;

% reshape the input field
n_vars = length(full_sol.variables);
n_nodes = size(reshape(full_sol.values(:,1), [], n_vars), 1);
p_sel = zeros(n_nodes, n_vars);
a_sel = ismember(full_sol.variables, {'a'});
p_sel(:, a_sel) = 1; p_sel = logical(reshape(p_sel, [], 1));

% assume slight compressibility of the fluid and convert pressure into
% temperature using a isentropic relation
temp_rel = @(a) a.^2 / (kappa*R) ;
is_rel = @(t) p_ref * (t/t_ref).^(kappa/(kappa-1)) - p_ref;

% decompose velocity and pressure field seperately
t_field = temp_rel(full_sol.values(p_sel, :));
full_sol.values(p_sel, :) = is_rel(t_field);

if any(a_sel)
    full_sol.variables{a_sel} = 'p'; 
else
   p_sel = zeros(n_nodes, n_vars);
   p_sel(:, ismember(full_sol.variables, {'p'})) = 1;
   p_sel = logical(reshape(p_sel, [], 1)); 
end
end

