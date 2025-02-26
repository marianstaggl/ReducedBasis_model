function [N, Nres] = get_applied_bc(self, N, Nres)
%GET_APPLIED_BC apply the boundary conditions at the problem
%   dirichlet boundaries: the equations are removed from the system as they
%   are not solved
%
%   neumann boundaries: these boundary conditions appear naturally on the
%   left-hand side of the problem. The actual value is imposed on the
%   righthandside
%
%   robin boundaries: these boundaries represent a linear combination of
%   the former two. an example is the convective boundary condition for the
%   heat transfer equation k*dq = h*(Tinf - Ti)

% robin: assume k and h is 1 and move the value to the left
r_sel = logical(self.rob_bounds(:,1));
N(r_sel, r_sel) = N(r_sel, r_sel) + speye(sum(r_sel));

% neumann: do nothing as they appear naturally

% dirichlet: remove the corresponding lines in the operator L
N = N(self.i_DoF, self.i_DoF); Nres = Nres(self.i_DoF);
end