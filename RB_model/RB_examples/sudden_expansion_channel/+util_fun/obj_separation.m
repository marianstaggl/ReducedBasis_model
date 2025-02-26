function bubble_length = obj_separation(ro_model, s)
% define some geometry parameters
n_query = 1000;
x_min = 5; x_max = 30;

% create the query line and evaluate the x velocity along this line
lower_pos = min(ro_model.fe_model.mesh.lines(1).pos(:,2));
x_line = linspace(x_min, x_max, n_query);
y_line = (lower_pos + 0.1) * ones(1, n_query);
prim_matrix = reshape(s, [], 3);
u_vel = ro_model.fe_model.mesh.get_fun_interp(...
    prim_matrix(:,2), [x_line', y_line']);

% take the sign of the velocity
bubble_fun = sigmoid(-1e6 * u_vel) - u_vel;
bubble_length = (x_max-x_min)/n_query * sum(bubble_fun);
end