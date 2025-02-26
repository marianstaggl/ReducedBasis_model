function fe_model = get_fe_model(~, root_folder)
%GET_FE_MODEL Create an FEmodel based on the files in the root folder
%   the finite element model is necessary to construct the 3rd order
%   tensors for the reduced model. as governing equations, the isentropic
%   navier stokes equations are chosen.

% create the fe model using the finite element mesh
fe_model = FEModel('isa. navier stokes');
fe_model.sys_variables = {'a', 'u', 'v'};
fe_model.set_mesh(fullfile(root_folder, 'mesh.m'), 1); % set the element
end

