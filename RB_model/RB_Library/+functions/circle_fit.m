function [x_c, y_c, radius, a] = circle_fit(x, y)
% fit a circle into a few points

% get the number of points
numPoints = numel(x);

% get the needed terms
xx = x .* x; yy = y .* y; xy = x .* y;

% form the matrix
A = [sum(x),  sum(y),  numPoints;
    sum(xy), sum(yy), sum(y);
    sum(xx), sum(xy), sum(x)];
B = [-sum(xx + yy) ;
    -sum(xx .* y + yy .* y);
    -sum(xx .* x + xy .* y)];

% calc the coefficients and the centers
a = A \ B; x_c = -.5 * a(1); y_c = -.5 * a(2);

% get the radius of the circle
radius  =  sqrt((a(1) ^ 2 + a(2) ^ 2) / 4 - a(3));

end