function u_plus = wake_law(y_plus, tau_w, u_mean)
%WALL_LAW produce a velocity profile according to the wake law
% define the parameters for wake law (channel-flow)
 rho = 1; u_tau = sqrt(abs(tau_w)/rho);
 d = 1e-2; d_plus = d*u_tau/14.6e-6;
 mue = 0.75*d_plus; sig = 0.25*d_plus;
 n = length(y_plus);

% generate the wall-law region
up_log = wall_law(y_plus);
up_fst = u_mean/u_tau*ones(1,n);

% blend the two profiles according to wake-law
ef = (1/2)*(1+erf((y_plus-mue)/sqrt(2*sig^2)));
u_plus = (1-ef).*up_log + ef.*up_fst;

end


function u_plus = wall_law(y_plus)
%WALL_LAW produce a velocity profile according to the wall law
% constant values from wikipedia
kappa = 0.384; c_plus = 4.17;

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
u_plus = makima(y_val, u_val, y_plus);
end