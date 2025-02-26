function new_res = interp_gappyPOD(bc_modes, query_prof, dir_sel)
%INTERP_GAPPYPOD use a gappyPOD approach to reconstruct the data
c_mean = bc_modes(:, 1); c_bound = bc_modes(:, 2:end); 

% use a ridge regression for the projection
X = c_bound(dir_sel, :);
y = (query_prof - c_mean(dir_sel));


c_proj = (X'*X + 1e-4*eye(size(c_bound, 2)))\(X'*y);
% c_proj = X\y;
new_res = c_mean + c_bound * c_proj; 
end
