function vel_prof = obj_velprofile(ro_model, s)
%OBJ_VELPROFILE match the velocity profile at a certain x position
% get the minimum and maximum values for the y positions
x_pos = 10;
y_min = min(ro_model.fe_model.mesh.nodes(:,2));
y_max = max(ro_model.fe_model.mesh.nodes(:,2));

% calculate the velocity amplitude of the vector s
prim_matrix = reshape(s, [], 3);
u_vel = prim_matrix(:,2);

% interpolate on a line from y_min to y_max
x_line = x_pos * ones(1, 200);
y_line = linspace(y_min, y_max, 200);
vel_prof = ro_model.fe_model.mesh.get_fun_interp(u_vel,...
    [x_line', y_line']);
end
