function b = get_b(self, nodes)
%GET_B Get the base function values at given points
%   this function returns the basefunction values of an element on the
%   given coordinates eta and zeta. It uses lagrange polynomials as
%   basefunctions.

% split up the points in eta and zeta coordinates
eta = nodes(:,1); zeta = nodes(:,2);

% get the supportpoints of the reference element
lag_support = linspace(-1,1,self.order+1);

% predefine the b array
b = zeros((self.order+1).^2, length(eta));

% loop over the support points and get the basefunctions
for i=1:(self.order+1)
    for j=1:(self.order+1)
        
        % get the lagrange polynomials
        Li = self.get_lagrange_base(eta, lag_support, i, 0);
        Lj = self.get_lagrange_base(zeta, lag_support, j, 0);
        
        % convert the subindices to basefunction index
        b_idx = sub2ind([self.order+1, self.order+1], i, j);
        
        % save the polynomial into the array
        b(b_idx,:) = reshape(Li.*Lj,1,[]);
    end
end
end