function [Re, Ma] = mean2nondim(u_m, t_m, L, p)
%MEAN2NONDIM Calculate the Reynolds and Machnumber

% get the machnumber
Ma = u_m/sqrt(1.4*287*t_m);

% define the constants for sutherlands law
mue_0 = 18.27e-6; T_0 = 291.15; C = 120;
mue = mue_0*(T_0 + C)/(t_m + C)*(t_m/T_0)^(3/2);

% get the reynoldsnumber
Re = (p*u_m*L)/(287*t_m*mue);
end

