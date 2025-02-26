function p = get_deim_points(U)
% This function carries out DEIM point selection using the columns of the matrix U.

% get the size of U
[~, n] = size(U);
p = zeros(n,1);

% Select the first DEIM point
[~, p(1)] = max(abs(U(:,1)));

% Identify additional DEIM points, reading in additional columns of U one at a time
for j = 2:n
    u = U(:,j);
    r = u - U(:,1:j-1)*(U(p(1:j-1),1:j-1)\u(p(1:j-1)));
    [~,p(j)] = max(abs(r));   
end
end