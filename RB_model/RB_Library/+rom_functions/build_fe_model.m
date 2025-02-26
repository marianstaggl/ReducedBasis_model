function fe_model = build_fe_model(eqn_type)
%BUILD_FE_MODEL setup the finite element model

if isequal(eqn_type, 'incompressible')
    fe_model = FEModel('inc. navier stokes'); % use the stokes equations
    fe_model.set_mesh('mesh.m', 1); % set the element
    fe_model.set_boundarys('boundary_inc.txt')

elseif isequal(eqn_type, 'isentropic')
    fe_model = FEModel('isa. navier stokes'); % use the stokes equations
    fe_model.set_mesh('mesh.m', 1); % set the element
    fe_model.set_boundarys('boundary_isa.txt')
end
end

