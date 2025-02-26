function [u_m,t_m] = nondim2mean(Re, Ma, L, p)
%NONDIM2MEAN use Re, Ma L and p to get the mean velocity and the mean
%temperature

% define the constants for sutherlands law
mue_0 = 18.27e-6; T_0 = 291.15; C = 120;

% calc t_m using Re and Ma
a = Re*mue_0*(T_0 + C)*287;
b = -Ma*sqrt(1.4*287)*T_0^(3/2)*p*L;
c = b*C;

t_m = roots([a,b,c]); t_m = t_m(t_m>0);

% calculate the u_m using sutherlands law
mue = mue_0*(T_0+C)/(t_m+C)*(t_m/T_0)^(3/2);
rho = p/(287*t_m);

u_m = (Re*mue)/(rho*L);
end

