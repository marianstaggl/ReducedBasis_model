function vel = apply_bl(vel, t_lo, d_lo, t_up, d_up)
%APPLY_BL apply a boundary layer to a vel profile
%   the velocity profile and the bl thickness is given to generate a bl

% reshape the vel profile
vel = reshape(vel, 1, []);

% define the length and velocity scale
l_sc = 0.1; vel_sc = 150;

% define the y vector
n = length(vel); y = linspace(0,l_sc,n);

% define the selfsimilar profiles at hub ans shroud
u_lo = vel(1)*((wake_law(y, t_lo, d_lo, vel_sc)/vel_sc)-1);
u_up = vel(end)*((wake_law(y, t_up, d_up, vel_sc)/vel_sc)-1);

% construct the full profile
if rem(length(y),2) ~= 0
    n = floor(length(y)/2);
    vel = vel + [u_lo(1:n+1), fliplr(u_up(1:n))];
else
    n = length(y)/2;
    vel = vel + [u_lo(1:n), fliplr(u_up(1:n))];
end
end

function u = wake_law(y, tau_w, d, u_mean)
%WALL_LAW produce a velocity profile according to the wake law
% define the parameters for wake law (channel-flow)
rho = 1; u_tau = sqrt(abs(tau_w)/rho);

% get the y_plus values
y_plus = y*u_tau/14.6e-6; d_plus = d*u_tau/14.6e-6;
mue = 0.75*d_plus; sig = 0.25*d_plus;
up_fst = u_mean/u_tau*ones(1,length(y));

% generate the wall-law region
up_log = wall_law(y_plus);

% blend the two profiles according to wake-law
ef = (1/2)*(1+erf((y_plus-mue)/sqrt(2*sig^2)));
u_plus = (1-ef).*up_log + ef.*up_fst;

% calculate the u value
u = u_plus*u_tau;
end

function u_plus = wall_law(y_plus)
% setup the parameters (these imply a B=4.17)
kappa = 0.384; a = -10.3061;
al = (-1/kappa-a)/2; be = sqrt(-2*a*al-al^2); 
R = sqrt(al^2 + be^2);

% define the composite profile (Musker 1979)
u_plus = 1/kappa*log((y_plus - a)/-a) + R^2/(a*(4*al-a))*((4*al + a)*...
    log(-a/R*sqrt((y_plus-al).^2 + be^2)./(y_plus - a)) + al/be*...
    (4*al+5*a)*(atan((y_plus-al)/be) + atan(al/be)));
end