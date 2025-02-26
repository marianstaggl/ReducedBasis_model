function c_mesh = static_get_c_mesh(b_ids, n_wall)
%GET_C_MESH Create a quadratic control mesh
%   this mesh is wrapped around the default calculation mesh. if needed
%   this mesh can be deformed to match a new geometry. all the nodes of the
%   calculation mesh inside this mesh are deformed accordingly

% define the number of wall points
ids = 1:1:n_wall; c_mesh = FEMesh();

% the mesh itself is a rectangle with small spacing on the wall side
x_l = linspace(-1, 1, n_wall); y_l = -1*ones(1, n_wall); id_l = b_ids(1)*ones(1, n_wall);
x_u = linspace(-1, 1, n_wall); y_u =  1*ones(1, n_wall); id_u = b_ids(2)*ones(1, n_wall);

% define the nodes array
c_mesh.nodes = [x_l', y_l'; x_u', y_u'];
c_mesh.elems = [ids(1:end-1)', ids(1:end-1)'+1, ids(1:end-1)' + n_wall, ids(1:end-1)' + n_wall + 1];
c_mesh.edges = vertcat([ids(1:end-1)', ids(1:end-1)'+1, id_l(1:end-1)'],...
    [ids(1:end-1)' + n_wall, ids(1:end-1)' + n_wall + 1, id_u(1:end-1)']);

% add elements to the mesh
c_mesh.ref_elem = FElement(1);
c_mesh.par_elem = FElement(1);

% set the lines of the fem model
c_mesh.lines = c_mesh.get_lines( c_mesh.nodes, c_mesh.edges );
end

