function set_mesh(self, meshfile, order)
%SET_MESH Summary of this function goes here
%   Detailed explanation goes here

% read the meshfile and set the corresponding properties
self.mesh = FEMesh( meshfile, order );

% setthe mapping functions of the mesh
self.set_mapping();
end

