function new_res = interp_hyb(bc_modes, par_mat, res_mat, query_par, query_prof, dir_sel)
%INTERP_HYB use a hybrid method to interpolate the result
% calculate the boundary amplitudes by projection
bound_sol = rom_functions.interp_gappyPOD(bc_modes, query_prof, dir_sel);

% calculate the inner modes by interpolations
inner_sol = rom_functions.interp_invD(par_mat, res_mat, query_par, 3, 1);

% sum up the result from boundary projection and inner interpolations
new_res = bound_sol + inner_sol;
end

