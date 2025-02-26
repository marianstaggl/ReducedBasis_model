function [bound_modes, domain_modes, dir_sel] = calc_boundary_modes(dir_sel, full_modes, cutoff)
% remove the mean and get the centered solution
mean_sol = mean(full_modes, 2);
cent_sol = full_modes - mean_sol;

% decompose the result only with respect to dirichlet bounds
[~, s, v_dec] = svd(cent_sol(dir_sel, :));
u_dec = cent_sol * v_dec;
u_dec = u_dec ./ vecnorm(u_dec, 2, 1);
bound_modes = [mean_sol, u_dec(:, 1:cutoff)];

% collect the remaining modes and orthogonalize them
[domain_modes, ~, ~] = svd(u_dec(:, cutoff+1:end ), 'econ');
end