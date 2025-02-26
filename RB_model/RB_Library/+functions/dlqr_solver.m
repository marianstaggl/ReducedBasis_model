function x = dlqr_solver(A, y)
% perform a qr decomposition
[q, r] = custom_qr_decomposition(A);

% solve the system by substitution
x = solve_upper_triangular(r, q'*y);
end

function [Q, R] = custom_qr_decomposition(A)
% preallocate some arrays
[m, n] = size(A);
R = A; Q = eye(m);

% check if A is a dlarray
if isa(A, 'dlarray'), Q = dlarray(Q); end

for k = 1:min(m, n)
    % Extract the column to work on
    x = R(k:end, k);

    % Compute the Householder vector
    v = x;
    v(1) = sign(x(1)) * sqrt(sum(x.^2)) + x(1); % sign(x(1)) * norm(x) + x(1);
    v = v / sqrt(sum(v.^2)); % norm(v);

    % Apply the transformation to R (vectorized)
    R(k:end, k:end) = R(k:end, k:end) - 2 * (v * (v' * R(k:end, k:end)));

    % Apply the transformation to Q (vectorized)
    Q(:, k:end) = Q(:, k:end) - 2 * (Q(:, k:end) * v) * v';
end

% return only the upper trianuglar matrix
R = R(1:n, :);
Q = Q(:, 1:n);
end

function x = solve_upper_triangular(U, b)
% Helper function to solve Ux = b for upper triangular U
n = size(U, 1);
x = zeros(n, 1);

% check if A is a dlarray
if isa(U, 'dlarray') | isa(b, 'dlarray'), x = dlarray(x); end

for i = n:-1:1
    x(i) = (b(i) - U(i, i+1:end) * x(i+1:end)) / U(i, i);
end
end