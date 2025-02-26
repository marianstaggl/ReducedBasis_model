function tri = get_triangulation(self)
%GET_TRIANGULATION Get a triangulated version of the mesh
%   Translate the quad mesh into a triangulated surface

% preallocate the arrays
ref_subs = self.ref_elem.get_sub_quads();
num_subs = size(ref_subs,1);
subs = zeros(self.n_elem * (self.ref_elem.order-1)^2, 4);

% loop over all of the elements
s_c = 0;
for i=1:self.n_elem
    % insert the subquads
    subs(s_c + 1:s_c + num_subs,:) = reshape(self.elems(i,ref_subs),[],4);
    s_c = s_c + num_subs;
end

% split the sub quads into triangles
tri = [subs(:, [1, 2, 3]); subs(:, [1, 3, 4])];
end

