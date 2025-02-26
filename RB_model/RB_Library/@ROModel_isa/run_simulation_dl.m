function curr_state = run_simulation_dl(self, prim_state, visc_state, ro_opts, n_iter)
%RUN_SIMULATION Summary run a simulation for the given boundary conditions
%   The equations of the reduced model are solved for the given boundary
%   conditions. number of iterations, degree of overdetermination, ... are
%   specified within the options struct.

% take the 4th order tensors and condense them into 3rd order
if isequal(ro_opts.build, 'approx'), self.set_csys_approx();
elseif isequal(ro_opts.build, 'exact'), self.set_csys_exact();
elseif isequal(ro_opts.build, 'keep') % do nothing in this case
else, error('the defined build option is not available...')
end

% get the current system matrices
diff_3rd = self.diff_sys.diff_3rd;
conv1_3rd = self.conv_sys.conv1_3rd;
conv2_3rd = self.conv_sys.conv2_3rd;

% specify the degrees of freedom of the system
curr_state = dlarray(prim_state);
[~, s_DoF, ~, t_DoF] = self.get_rom_base(ro_opts);

% solve the non linear system iteratively
for i=1:n_iter
    % get the systems equations
    [N, R] = get_linear_sys(diff_3rd, visc_state,...
        conv1_3rd, conv2_3rd, curr_state);
    N = N(t_DoF, s_DoF); R = R(t_DoF);

    % take the next step with underrelaxation
    dR = - 0.5*functions.dlqr_solver(N, R);
    
    % take a picard step
    curr_state(s_DoF) = curr_state(s_DoF) + dR;
end
end

function [Nr, Rr] = get_linear_sys(diff_3rd, curr_visc,...
    conv1_3rd, conv2_3rd, curr_state)
%GET_LINEAR_SYS get the linearized system around the current state
%   the current state is used to assemble the linearized system around it.
%   To achieve this, the 3rd order tensors are contracted along the k index.
Nr_visc = sum(diff_3rd.*reshape(curr_visc,1,1,[]),3);
Nr1_conv = sum(conv1_3rd.*reshape(curr_state,1,1,[]),3);
Nr2_conv = sum(conv2_3rd.*reshape(curr_state,1,1,[]),3);

% build the reduced system
Nr1 = Nr1_conv + Nr_visc;
Rr = Nr1 * curr_state;
Nr = Nr1 + Nr2_conv; % add the additional terms to get the full newton method
end

function plot_residuals(res,counter)
% plot the current residuals
clf; semilogy(res); text(1, 1, ['curr res: ' num2str(res(counter)) '  '],...
    'Units','normalized', 'HorizontalAlignment','right',...
    'VerticalAlignment','top');

% grid on and draw now
grid on; drawnow;
end