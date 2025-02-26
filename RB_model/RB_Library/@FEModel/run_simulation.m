function [curr_state, res, dres, time_record] = run_simulation(self, n_iter, opts)
%SOLVE_SYSTEM Solve the state equations in an iterative manner
%   there are various options to run a simulation. normally for a navier
%   stokes equation one would run a few picard iteration steps (setting the
%   option 'lin' to 0 and the converge the system with a newtons method by
%   setting 'lin' to 1 again.

% get the initial conditions
sel = logical(self.dir_bounds(:,1));
curr_state = self.get_init_cond(opts);

% get the righthandside vector
f = self.get_rhs_vec();
curr_state(sel) = f(sel);

% if a record should be made
if opts.time_rec, time_record = zeros(length(curr_state), n_iter); 
else, time_record = nan; end

% run the calulation in a loop
for i=1:n_iter
    
    % get the systems equations
    [N, R, c_res, curr_state] = self.get_linear_sys(curr_state, opts,...
        self.sys_equations);
    
    % apply the boundary conditions
    [N, R] = self.get_applied_bc(N, (R-f));
    
    % take the next step
    dR = - opts.damping*(N\R);
    
    % take a picard step
    curr_state(self.i_DoF) = curr_state(self.i_DoF) + dR;
    
    % save the current residual and the stepsize
    res(i) = c_res; dres(i) = norm(dR);
    
    % store the current step if necessary
    if opts.time_rec, time_record(:,i) = curr_state; end
    
    % plot the residuals
    if opts.plot, plot_residuals(res, i); end
end
% return only the last residual
res = res(end); dres = dres(end);

% plot the result
if opts.plot, self.plot_field(curr_state); end
end

function plot_residuals(res,counter)

% clear the figure
clf;

% plot the current residuals
semilogy(res); text(1, 1, ['curr res: ' num2str(res(counter)) '  '],...
    'Units','normalized', 'HorizontalAlignment','right',...
    'VerticalAlignment','top');

% grid on and draw now
grid on; drawnow;
end

