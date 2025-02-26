function [x_vel, y_vel] = calc_blasius_solution(x_vec, y_vec, u_inf, nue)
% Use fsolve to ensure the boundary function is zero. The result is the
% unknown initial condition.
opt = optimset('Display','off','TolFun',1E-20);
F = fsolve(@(F) eval_boundary(F),[0,0,0.33],opt);

% Solve the ODE-IVP with the converged initial condition
[x,y] = solve_ode(F);

% use the blasius solution to calculate the true x and y velocities
[X, Y] = meshgrid(x_vec, y_vec);
eta = Y .* sqrt(u_inf ./ (nue * X));

% evaluate the blasius solution at the eta points
f_eta = interp1(x, y(:,1), eta, 'linear', 'extrap');
fs_eta = interp1(x, y(:,2), eta, 'linear', 'extrap');

% calculate y velocity
x_vel = u_inf * fs_eta;
y_vel = 0.5 * sqrt((nue*u_inf)./X) .* (eta .* fs_eta - f_eta);
end

function [x,y] = solve_ode(F)
% Solve the ODE-IVP with initial condition F on [0 100] (arbitrary upper
% bound)
[x,y] = ode45(@(x,y) [y(2); y(3); -0.5*y(1)*y(3)],[0 300],F); %solve BVP                
end

function [g] = eval_boundary(F)
% Get the solution to the ODE with inital condition F
[x,y] = solve_ode(F);

% Get the function values (for BCs) at the starting/end points
f_start = y(1,1); %f(0) = 0
df_start = y(1,2); %f'(0) = 0
df_end = y(end,2); %f'(inf) - 1 = 0

% Evaluate the boundary function
g = [f_start
     df_start
     df_end - 1];
end