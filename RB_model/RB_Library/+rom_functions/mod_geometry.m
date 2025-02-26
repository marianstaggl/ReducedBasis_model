function fe_model = mod_geometry(c_nodes, fe_model, par_vec)
% modify the current node positions
new_nodes = c_nodes; x_vec = new_nodes(:, 1) - min(new_nodes(:, 1));
new_nodes(:,2) = new_nodes(:,2) + par_vec(1)*x_vec.^2 + par_vec(2)*x_vec.^3;

% change them in the model itself
fe_model.mesh.set_nodes(new_nodes); 
fe_model.set_mapping();
end