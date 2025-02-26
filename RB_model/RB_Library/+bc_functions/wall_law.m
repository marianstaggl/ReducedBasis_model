function [y_val,u_val] = wall_law(n)
%WALL_LAW produce a velocity profile according to the wall law
% constant values from wikipedia
kappa = 0.41; c_plus = 5;

% get the visc law and it's domain
visc_dom = linspace(0, 5, 1e4);
visc_law = visc_dom;

% get the log law region and its domain
log_dom = linspace(30, 700, 1e4);
log_law = 1/kappa*log(log_dom) + c_plus; 

% get the outer domain with constant value
out_dom = linspace(950, 1000, 1e4);
out_law = (1/kappa*log(850) + c_plus) * ones(1,length(out_dom));

% get the piecewise profile
y_plus = [visc_dom, log_dom, out_dom];
u_plus = [visc_law, log_law, out_law];

% use a makima spline to interpolate the profile
y_val = linspace(0, 1000, n);
u_val = makima(y_plus, u_plus, y_val);
end

