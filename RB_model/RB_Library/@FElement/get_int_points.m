function [points, weights] = get_int_points( order )
%GET_INT_POINTS Get the integration points and their weights
%   for the gauss quadrature the position of the integration points as well
%   as their weights are needed. the summed weights (transformed if needed)
%   give the area of the element

% get the onedimensional gausian points and weights
[xi, wi] = gaussquadrules(order+1);

% scale up to second dimension
[ZETA, ETA] = meshgrid(xi, xi);
[WETA, WZETA] = meshgrid(wi, wi);

% reshape the 2 arrays
points = [reshape(ETA,[],1), reshape(ZETA,[],1)];
weights = reshape(WETA.*WZETA, 1, []);

end

function [x, w] = gaussquadrules(n)
% Generates the abscissa and weights for a Gauss-Legendre quadrature.
% Reference:  Numerical Recipes in Fortran 77, Cornell press.
x = zeros(n,1);                                           % Preallocations.
w = x;
m = (n+1)/2;
for ii=1:m
    z = cos(pi*(ii-.25)/(n+.5));                        % Initial estimate.
    z1 = z+1;
while abs(z-z1)>eps
    p1 = 1;
    p2 = 0;
    for jj = 1:n
        p3 = p2;
        p2 = p1;
        p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj;       % The Legendre polynomial.
    end
    pp = n*(z*p1-p2)/(z^2-1);                        % The L.P. derivative.
    z1 = z;
    z = z1-p1/pp;
end
    x(ii) = -z;                                   % Build up the abscissas.
    x(n+1-ii) = z;
    w(ii) = 2/((1-z^2)*(pp^2));                     % Build up the weights.
    w(n+1-ii) = w(ii);
end
end