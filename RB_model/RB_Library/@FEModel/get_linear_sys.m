function [N,R,res,curr_state] = get_linear_sys(self, curr_state, opts, sys_eqn)
%GET_LINEAR_SYS Get a linearized version of the systems equations
%   depending on the choosen system of equations a different function is
%   used to assemble the global stiffness matrix. up until now the heat
%   transfer equation, the stokes equation and the navier stokes equations
%   are available

% change between stokes/navier stokes/heat
switch sys_eqn
    case 'potential flow'
        % get the linear system of heat transfer
        [N, R] = self.get_OP_PotentialFlow(curr_state, opts);
        res = norm(R(self.i_DoF));
        
    case 'stokes equations'
        % set the density to 0 to switch off the convective terms
        [N, R] = self.get_OP_StokesFlow(curr_state, opts); res = norm(R);
        
    case 'inc. navier stokes'
        % get the galerkin part of the matrix
        [N1, N2, R] = self.get_OP_incNavierStokes(curr_state, opts);
        res = norm(R(self.i_DoF));

        % and create the linearized system
        N = (N1 + opts.lin*N2);

        % get the SUPG part of the matrix
        if opts.stab > 0
            [Ns1, Ns2, Rs] = self.get_OP_incNavierStokes_SUPG(curr_state, opts);

            % form the global matrix
            N = N + opts.stab*(Ns1 + opts.lin*Ns2);
            R = R + opts.stab*Rs;
        end
        
    case 'ish. navier stokes'
        % get the galerkin part of the matrix
        [N1, N2, R] = self.get_OP_ishNavierStokes(curr_state, opts);
        res = norm(R(self.i_DoF));

        % and create the linearized system
        N = (N1 + opts.lin*N2);

        % get the SUPG part of the matrix
        if opts.stab > 0
            [Ns, Rs] = self.get_OP_ishNavierStokes_SUPG(curr_state, opts);

            % form the global matrix
            N = N + opts.stab*Ns;
            R = R + opts.stab*Rs;
        end

    case 'isa. navier stokes'
        % get the galerkin part of the matrix
        [N1, N2, R] = self.get_OP_isaNavierStokes(curr_state, opts); 
        res = norm(R(self.i_DoF));

        % and create the linearized system
        N = (N1 + opts.lin*N2);

        % get the SUPG part of the matrix
        if opts.stab > 0
            [Ns, Rs] = self.get_OP_isaNavierStokes_SUPG(curr_state, opts);

            % form the global matrix
            N = N + opts.stab*Ns;
            R = R + opts.stab*Rs;
        end

    case 'turb. isa. navier stokes'
        % get the turbulence field for the current flow field
        [opts, curr_state] = self.get_Turbulence_field_SA(curr_state, opts);

        % get the galerkin part of the matrix
        [N1, N2, R] = self.get_OP_isaNavierStokes(curr_state, opts); 
        res = norm(R(self.i_DoF));
        
        % get the SUPG part of the matrix
        [Ns, Rs] = self.get_OP_isaNavierStokes_SUPG(curr_state, opts);
        
        % form the global matrix
        N = (N1 + opts.lin*N2) + opts.stab*Ns; R = R + opts.stab*Rs;
    
    case 'turb. inc. navier stokes'
        % get the turbulence field for the current flow field
        [opts, curr_state] = self.get_Turbulence_field_SA(curr_state, opts);

        % get the galerkin part of the matrix
        [N1, N2, R] = self.get_OP_incNavierStokes(curr_state, opts); 
        res = norm(R(self.i_DoF));
        
        % get the SUPG part of the matrix
        [Ns, Rs] = self.get_OP_incNavierStokes_SUPG(curr_state, opts);
        
        % form the global matrix
        N = (N1 + opts.lin*N2) + opts.stab*Ns; R = R + opts.stab*Rs;
        
    otherwise
        error('no valid equation was specified')
        
end
end

