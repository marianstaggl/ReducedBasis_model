function b = get_dbdz(self, nodes, d_order)
%GET_DBZ Get the base function zeta derivatives at given points
%   this function returns the basefunction zeta derivatives of the reference
%   element the given coordinates eta and zeta. It uses lagrange 
%   polynomials as basefunctions.

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
        Lj = self.get_lagrange_base(zeta, lag_support, j, d_order);
        
        % convert the subindices to basefunction index
        b_idx = sub2ind([self.order+1, self.order+1], i, j);
        
        % save the polynomial into the array
        b(b_idx,:) = reshape(Li.*Lj,1,[]);
    end
end
end