function sys = build_rom_system(fe_model, full_shape, full_test)
opts = fe_model.get_opts('lin', 0, 'stab', 0);

% preallocate array sizes
ndb = size(full_shape); ndt = size(full_test);
sys = zeros(ndt(2), ndb(2), ndb(2)+1);

% get the diffusive part of the system
[N_diff, ~, ~] = fe_model.get_linear_sys(zeros(ndb(1),1), opts, ...
    fe_model.sys_equations);
sys(:,:,1) = full_test' * (N_diff*full_shape);

% loop over the modes and get the convective parts
for i=1:ndb(2)
    % show the progress
    disp(['building matrix ', num2str(i), ' out of ', num2str(ndb(2))]);
    
    % use each mode to create the matrix
    [N_conv, ~, ~] = fe_model.get_linear_sys(full_shape(:,i), opts, ...
        fe_model.sys_equations);
    
    % project them onto the modes
    sys(:,:,i+1) = full_test'*(N_conv-N_diff)*full_shape;
end
end