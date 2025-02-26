function init_cond = get_init_cond(self, opts)
%GET_INIT_COND get the initial conditions for nonlinear PDEs
%   to solve nonlinear pdes like the navier stokes equation we have to make
%   an initial guess to start from (with iterative methods). the initial
%   guess for the navier stokes equations is calculated by a stokes
%   solution.

% check if there is a initial condition passed
if ~isnan(opts.init)
    init_cond = opts.init; return;
end

% see what the systems equations are and then initialize
switch self.sys_equations
    case 'potential flow'
        init_cond = init_heat_transfer(self, opts);
    case 'stokes equations'
        init_cond = init_stokes_eqn(self, opts);
    case 'inc. navier stokes'
        init_cond = init_inc_navier_stokes(self, opts);
    case 'turb. inc. navier stokes'
        init_cond = init_turb_inc_navier_stokes(self, opts);
    case 'isa. navier stokes'
        init_cond = init_isa_navier_stokes(self, opts);
    case 'turb. isa. navier stokes'
        init_cond = init_isa_navier_stokes_turb(self, opts);
    otherwise
        error('the specified equations are not available...')
end
end

function init_cond = init_isa_navier_stokes_turb(self, opts)
% don't solve the nt equations in the normal system
i_DoF = reshape(self.i_DoF, [], numel(self.sys_variables));
i_DoF(:, ismember(self.sys_variables, 'nt')) = false;
self.i_DoF = reshape(i_DoF, [], 1);

% get the current initial conditions and set the viscosity very high
init_cond = reshape(self.dir_bounds(:,2), [], numel(self.sys_variables));
dir_sel = reshape(self.dir_bounds(:,1), [], numel(self.sys_variables));
asel = ismember(self.sys_variables, 'a'); amax = max(init_cond(:,asel));
usel = ismember(self.sys_variables, 'u'); umax = max(init_cond(:,usel));

% get an incompressible and viscous flow regime
true_avec = init_cond(:, asel); true_avec(~dir_sel(:,asel)) = amax;
init_cond(:, asel) = 100*umax; opts.mue = 100*opts.mue;

% get the calculate the flow field
[L,R] = self.get_linear_sys(init_cond, opts, 'isa. navier stokes');
L = L(self.i_DoF,self.i_DoF); R = R(self.i_DoF);
init_cond(self.i_DoF) = init_cond(self.i_DoF) - L\R;

% return the initial conditions
init_cond(:, asel) = true_avec;
init_cond = reshape(init_cond, [], 1);
nt_sel = (~self.dir_bounds(:,1) & ~self.i_DoF);
init_cond(nt_sel) = 5*opts.mue/100;
end

function init_cond = init_isa_navier_stokes(self, opts)
% get the current initial conditions and set the viscosity very high
init_cond = reshape(self.dir_bounds(:,2), [], numel(self.sys_variables));
dir_sel = reshape(self.dir_bounds(:,1), [], numel(self.sys_variables));
asel = ismember(self.sys_variables, 'a'); amax = max(init_cond(:,asel));
usel = ismember(self.sys_variables, 'u'); umax = max(init_cond(:,usel));

% get an incompressible and viscous flow regime
true_avec = init_cond(:, asel); true_avec(~dir_sel(:,asel)) = amax;
init_cond(:, asel) = 100*umax; opts.mue = 100*opts.mue;

% get the calculate the flow field
[L,R] = self.get_linear_sys(init_cond, opts, 'isa. navier stokes');
L = L(self.i_DoF,self.i_DoF); R = R(self.i_DoF);
init_cond(self.i_DoF) = init_cond(self.i_DoF) - L\R;

% return the initial conditions
init_cond(:, asel) = true_avec;
init_cond = reshape(init_cond, [], 1);
end

function init_cond = init_turb_inc_navier_stokes(self, opts)
% initialize the turbulent flow with a stokes solution

% get the initial conditions
init_cond = self.dir_bounds(:,2);

% modify the inner degrees of freedom (don't solve the turb eqn.)
i_DoF = reshape(self.i_DoF, [], numel(self.sys_variables));
i_DoF(:, ismember(self.sys_variables, 'nt')) = false;
self.i_DoF = reshape(i_DoF, [], 1);

% get the stokes operator for the current system by setting rho to 0
[L,R] = self.get_linear_sys(init_cond, opts, 'stokes equations');

% calculate only the degrees of freedom (no dir bounds)
L = L(self.i_DoF, self.i_DoF); R = R(self.i_DoF);

% set the inner nodes as a solution to the stokes problem
init_cond(self.i_DoF) = init_cond(self.i_DoF) - L\R;
nt_sel = (~self.dir_bounds(:,1) & ~self.i_DoF);
init_cond(nt_sel) = 5*opts.mue; % set the freestream turbulence
end

function init_cond = init_inc_navier_stokes(self, opts)
% initialize navier stokes equations with a stokes solution

% get the current initial conditions
init_cond = self.dir_bounds(:,2);

% get the stokes operator for the current system by setting rho to 0
[L,R] = self.get_linear_sys(init_cond, opts, 'stokes equations');

% calculate only the degrees of freedom (no dir bounds)
L = L(self.i_DoF,self.i_DoF); R = R(self.i_DoF);

% set the inner nodes as a solution to the stokes problem
init_cond(self.i_DoF) = init_cond(self.i_DoF) - L\R;

end

function init_cond = init_stokes_eqn(self, ~)
% initial conditions of the stokes problem are simply the dirichlet bounds
init_cond = self.dir_bounds(:,2);
end

function init_cond = init_heat_transfer(self, ~)
% initial conditions of the heat transfer are simply the dirichlet bounds
init_cond = self.dir_bounds(:,2);
end

%{
function init_cond = init_isa_navier_stokes_turb(self, ~)
% initial conditions of the isentropic problem are simply the dirichlet bounds
init_cond = reshape(self.dir_bounds(:,2),[],4);
sel = reshape(self.dir_bounds(:,1),[],4);

% set the enthalpy and the u velocity everywhere
init_cond(:,1) = max(init_cond(:,1));
init_cond(~sel(:,2),2) = max(init_cond(:,2));
init_cond(~sel(:,3),3) = max(init_cond(:,3));
init_cond = reshape(init_cond,[],1);
end

function init_cond = init_isa_navier_stokes(self, ~)
% initial conditions of the isentropic problem are simply the dirichlet bounds
init_cond = reshape(self.dir_bounds(:,2),[],3);
sel = reshape(self.dir_bounds(:,1),[],3);

% set the enthalpy and the u velocity everywhere
init_cond(:,1) = max(init_cond(:,1)); init_cond(~sel(:,2),2) = max(init_cond(:,2));
init_cond = reshape(init_cond,[],1);
end
%}
