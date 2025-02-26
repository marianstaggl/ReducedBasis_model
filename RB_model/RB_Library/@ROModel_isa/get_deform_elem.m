function deim_elems = get_deform_elem(self, n_marker, deim_ids)
%GET_DEFORM_ELEM_ID get the element id of the deim points

% get the corresponding element for each DEIM point
num = size(self.fe_model.mesh.weights);
c_vec = reshape(repmat((1:1:num(1))', num(2)), [], 1);
elem_ids = reshape(repmat(c_vec, n_marker), [], 1);
deim_elems = elem_ids(deim_ids);
end

