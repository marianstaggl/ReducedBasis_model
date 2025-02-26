function set_csys_exact(self)
% preallocate array sizes
ndb = size(self.prim_shape); 
ndt = size(self.prim_test); 
ndf = size(self.visc_shape);

% preallocate the different tensors
conv1_sys = zeros(ndt(2), ndb(2), ndb(2)); % tensors for the convective parts
conv2_sys = zeros(ndt(2), ndb(2), ndb(2)); % tensors for the convective parts
diff_sys = zeros(ndt(2), ndb(2), ndf(2)); % tensors for the diffusive parts

% get the diffusive part of the 0 system
physics_opts = self.fe_model.get_opts('mue', 0);

% loop over the modes and get the convective parts
for i=1:ndb(2)
    % show the progress
    disp(['building convective matrix ', num2str(i), ' out of ', num2str(ndb(2))]);
    
    % use each mode to create the matrix
    [N1_conv, N2_conv, ~] = self.fe_model.get_OP_isaNavierStokes_conv(...
        self.prim_shape(:,i), physics_opts);
    
    % project them onto the modes
    conv1_sys(:,:,i) = self.prim_test'*N1_conv*self.prim_shape;
    conv2_sys(:,:,i) = self.prim_test'*N2_conv*self.prim_shape;
end

% loop over the turbulence fields and get the diffusive parts
for i=1:ndf(2)
    % show the progress
    disp(['building diffusive matrix ', num2str(i), ' out of ', num2str(ndf(2))]);

    % change the turbulence field of the current stage
    physics_opts.mtvec = self.visc_shape(:, i);
    [N_diff ,~] = self.fe_model.get_OP_isaNavierStokes_diff(...
        zeros(ndb(1),1), physics_opts);

    % project them onto the modes
    diff_sys(:,:,i) = self.prim_test'*N_diff*self.prim_shape;
end

% set the current system matrices
self.diff_sys.diff_3rd = diff_sys;
self.conv_sys.conv1_3rd = conv1_sys;
self.conv_sys.conv2_3rd = conv2_sys;
end