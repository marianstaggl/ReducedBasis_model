function rhs_vector = get_rhs_vec(self)
%GET_RHS_VEC Returns the right hand side vector of the current system
%   the rhs vectors "contains" the boundary conditions of the current
%   system and allows for the solution of the linearized version of the
%   system.
%
%   NOTE!! up til now, the neumann bounds do NOT get integrated even if
%   they should be. the reason is, that for the current aim all of the
%   neumann bounds will have a value of zero resulting also in an integral
%   of 0! If someone want's to use nonzero neumann bounds one has to
%   perform an integration if the grid on the bound is not equally spaced!
%   
%   NOTE!! the current conditions at an outflow of the domain are so called
%   do nothing boundary condition. we just set the forcing term to zero and
%   do nothing beside this.

% get the dirichlet and neumann boundary conditions
rhs_vector = self.dir_bounds(:,2) + self.neu_bounds(:,2) + self.rob_bounds(:,2);

end

