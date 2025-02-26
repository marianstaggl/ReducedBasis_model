function [y_val,u_val] = power_law(n, n_exp)
%POWER_LAW function to create a velocity profile accoring to a power law

% create a y vector ranging from 0 to 1
y_val = linspace(0,1,n); s_vec = interp1([0, 0.5, 1], [0, 0, 1], y_val)';

% create the velocity profile according to the power law
for i=1:length(n_exp), u_val(:,i) = (1-s_vec).^(1./n_exp(i)); end

end

