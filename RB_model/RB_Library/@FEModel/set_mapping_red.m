function set_mapping_red( self, elem_ids )
%SET_MAPPING set the mapped derivatives from the mesh object
%   to define the weak form of the state equations one needs to have the
%   derivatives of the base functions at each quadrature point. to get them
%   we map the basefunctions defined at the reference element onto the real
%   elements.

% set the weights and derivatives
[self.weights, self.dbdx, self.dbdy, self.d2bdx2, self.d2bdy2,...
    self.inv_jac, self.dbdx_ref, self.dbdy_ref] = self.mesh.get_map_dl(elem_ids);

% set the base function
ref_elem = self.mesh.ref_elem;
self.b = ref_elem.get_b(ref_elem.int_points);
end

