function Lj = get_lagrange_base(eta, eta_support, base_idx, derivation_order)
%GET_LAGRANGE_BASE function to calculate lagrange base polynomials
%   the lagrange ith polynomial are calculated on the given eta vector and
%   the order is given by the length of the eta support wich defines the
%   support points for the polynomial. If needed also derivatives of the
%   functions are assembled.

% enforce columnvec for eta and rowvec for eta_support to match dims
eta = reshape(eta,[],1); eta_support = reshape(eta_support,1,[]);

% get the single components of the function
Lj = (eta - eta_support)./(eta_support(base_idx) - eta_support);
dj = 1./(eta_support(base_idx) - eta_support);

% get the product of the single components
rem_idx = setdiff(1:1:size(dj,2), base_idx);
Lj = derive(Lj, dj, rem_idx, derivation_order)';

end

function dL = derive(Lj, dj, rem_idx, order)
% use a recursive function to get the derivative of the lagrange polynomial
% when the final stage is reached, the polynomial itself is returned,
% otherwise we call the function again

% check the order of the desired derivation
if (order ~= 0)
    % preallocate the array
    dL = zeros(size(Lj,1),1);
    
    % loop to get the sum of the components
    for l=1:length(rem_idx)
        
        % go one step deeper
        new_idx = setdiff(rem_idx, rem_idx(l));
        dL = dL + dj(:,rem_idx(l)).*derive(Lj, dj, new_idx, order-1);
    end
    
else
    % if the derivation order is 0 the derivative is the function itself
    dL = prod(Lj(:, rem_idx),2);
end

end