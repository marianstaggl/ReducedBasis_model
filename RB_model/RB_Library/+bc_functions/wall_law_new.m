function [y_val,u_val] = wall_law_new(n, tau_w, u_mean)
%WALL_LAW produce a velocity profile according to the wall law
% constant values from wikipedia
kappa = 0.41; c_plus = 5; 
m_vl = 15; sig_vl = 10;
m_lf = 5e3; sig_lf = 3e3;
n = 1e3;

% define the y+ values and the interfaces
y_plus = linspace(0, 1e4, 1e4);
E_vl = (1/2)*(1+erf((y_plus-m_vl)/sqrt(2*sig_vl^2))); % blending for laminar and log law
E_lf = (1/2)*(1+erf((y_plus-m_lf)/sqrt(2*sig_lf^2))); % blending for log law and freestream (aka law of the wake)

% calculate the u+ values of the profile
visc_law = y_plus; % viscosity in the sublayer
log_law = 1/kappa*log(y_plus) + c_plus; % velocity due to log law
free_stream = (u_mean/tau_w)*ones(1,n); % dimensionless freestream velocity

% define the self-similar velocity profile
u_plus = visc_law*(1-E_vl) +... % using the visc sublayer
    E_vl*log_law*(1-E_lf) +... % the log-law region
    E_lf*free_stream; % and the freestream region

% get the piecewise profile
y_plus = [visc_dom, log_dom, out_dom];
u_plus = [visc_law, log_law, out_law];

% use a makima spline to interpolate the profile
y_val = linspace(0, 1000, n);
u_val = makima(y_plus, u_plus, y_val);
end

function [y_plus,u_plus] = wall_law(n)
%WALL_LAW produce a velocity profile according to the wall law
% constant values from wikipedia
kappa = 0.41; c_plus = 5;

% define the domains of visc and log law
visc_dom = linspace(0, 5, 5); % domain of the visc-sublayer
log_dom = linspace(30, 1e4, 1e3); % reaches into the freestream

% define the u_plus values of  visc and log
visc_law = visc_dom;
log_law = 1/kappa*log(log_dom) + c_plus;

% get the piecewise profile
y_val = [visc_dom, log_dom];
u_val = [visc_law, log_law];

% interpolate with makime for buffer layer
y_plus = linspace(0, 1e4, n);
u_plus = makima(y_val, u_val, y_plus);
end

