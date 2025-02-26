function par2textfile(exp_folder, par, r_vec)
%PAR2TEXTFILE export a parametrization as a textfile

% sort the r_vector in ascending order and remove double values
r_vec = reshape(unique(sort(r_vec)),[],1); 

% get a normalized version of the r vector
r_n = ((r_vec - min(r_vec))/(max(r_vec) - min(r_vec)));

% get the different profiles for c_ax, c_circ and c_rad
[c_ax, c_circ, c_rad] = bc_functions.par2profile(par, 1e4);

% resample all the profiles on the given r vector
c_ax = interp1(linspace(0,1,1e4), c_ax, r_n);
c_circ = interp1(linspace(0,1,1e4), c_circ, r_n);
c_rad = interp1(linspace(0,1,1e4), c_rad, r_n);

% get the data into a table
T = table(r_vec, c_ax, c_circ, c_rad);
writetable(T, fullfile(exp_folder, 'velocity_profiles.csv'));

% make a file for the total temperature and the static pressure
dlmwrite(fullfile(exp_folder, 'total_temperature_in.txt'), par.t_tot_in);
dlmwrite(fullfile(exp_folder, 'static_pressure_out.txt'), par.ps_out);
end

