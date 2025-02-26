function [diff_sys, conv_sys] = get_sys_build(self, opts)
%SET_SYS_MAT create the 3rd order tensors for the reduced system
%   The finite element system is used to set up the 3rd order tensors for
%   the reduced model. One tensor is created for the diffusive part and
%   another tensor is created for the convective part. The 4th index of the
%   tensors is reserved for geometrical variations

if isempty(opts.geo_variations), geo_variations = {self.def_geo.def_nodes};
else, geo_variations = opts.geo_variations; end

% preallocate the necessary cell arrays
n = numel(geo_variations);
all_elem = 1:1:self.fe_model.mesh.n_elem;
diff_sys_cell = cell(1, n); 
conv1_sys_cell = cell(1, n);
conv2_sys_cell = cell(1, n);

deform_diff = cell(1, n); 
deform_conv = cell(1, n);

% loop over the different geometries and build multiple systems
node_size = size(self.fe_model.mesh.nodes);
for i=1:numel(geo_variations)
    % extract the nodes from the cell array
    new_nodes = geo_variations{i};
    if size(new_nodes) ~= node_size
        error(['the ' num2str(i) 'th geometry does not match in size\n' ...
            'the given size is ' num2str(size(new_nodes)) 'instead of ' ...
            num2str(node_size)])
    end
    
    % deform the geometry and get the system tensors
    self.fe_model.mesh.set_nodes(new_nodes); 
    self.fe_model.set_mapping(); self.set_csys_exact();
    diff_sys_cell{i} = self.diff_sys.diff_3rd;
    conv1_sys_cell{i} = self.conv_sys.conv1_3rd;
    conv2_sys_cell{i} = self.conv_sys.conv2_3rd;
    
    % get the deformation functions for the different geometries
    [deform_diff{i}, diff_marker] = self.get_deform_diff(all_elem);
    [deform_conv{i}, conv_marker] = self.get_deform_conv(all_elem);
end

% convert the cell arrays to 4th order tensors
[diff_sys_q, diff_xi, diff_deim] = build_system(...
    full(cell2mat(deform_diff)), {cat(4, diff_sys_cell{:})},...
    opts.n_geo_diff, 10 * opts.n_geo_diff, opts.plot);
[conv_sys_q, conv_xi, conv_deim] = build_system(...
    full(cell2mat(deform_conv)), {cat(4, conv1_sys_cell{:}),...
    cat(4, conv2_sys_cell{:})}, opts.n_geo_conv,...
    10 * opts.n_geo_conv, opts.plot);

% pack the results into structs and return them
diff_sys.tensor = diff_sys_q{1};
diff_sys.diff_xi = diff_xi;
diff_sys.diff_deim = diff_deim;
diff_sys.deim_elem = self.get_deform_elem(diff_marker, diff_deim);

conv_sys.tensor_1 = conv_sys_q{1};
conv_sys.tensor_2 = conv_sys_q{2};
conv_sys.conv_xi = conv_xi;
conv_sys.conv_deim = conv_deim;
conv_sys.deim_elem = self.get_deform_elem(conv_marker, conv_deim);
end

function [A_q, U_masked, deim_pts] = build_system(xi_mat, sys_mat,...
    n_base, n_deim, plot_switch)
% perform svd on the deformation matrices
xi_mat(any(isnan(xi_mat')), :) = 0;
[U, S, V] = svd(xi_mat, 'econ');
U = U .* diag(S)';

% plot the decay of the singular values
if plot_switch
    semilogy(diag(S(1:end-1, 1:end-1)))
    xlabel('number of bases')
    ylabel('singular value')
    grid on
end

% identify the DEIM points of the matrix
n_deim = min([n_deim, size(U, 2)]);
n_base = min([n_base, size(U, 2)]);
deim_pts = rom_functions.get_deim_points(U); 
deim_pts = deim_pts(1:n_deim);
U_masked = U(deim_pts, 1:n_base);

% create the 3rd order tensors
sa = size(sys_mat{1}); sa(4) = n_base; sa(sa==0) = 1; 
A_q = repmat({zeros(sa)}, 1, numel(sys_mat));
for i = 1:n_base
    for j = 1:numel(sys_mat)
        A_q{j}(:, :, :, i) = sum(sys_mat{j}.*reshape(V(:, i),1,1,1,[]),4);
    end
end
end

function [diff_sys, conv1_sys, conv2_sys] = get_single_geo_sys(c_rom)
% preallocate array sizes
ndb = size(c_rom.prim_shape); 
ndt = size(c_rom.prim_test); 
ndf = size(c_rom.visc_shape);

% preallocate the different tensors
conv1_sys = zeros(ndt(2), ndb(2), ndb(2)); % tensors for the convective parts
conv2_sys = zeros(ndt(2), ndb(2), ndb(2)); % tensors for the convective parts
diff_sys = zeros(ndt(2), ndb(2), ndf(2)); % tensors for the diffusive parts

% get the diffusive part of the 0 system
physics_opts = c_rom.fe_model.get_opts('mue', 0);

% loop over the modes and get the convective parts
for i=1:ndb(2)
    % show the progress
    disp(['building convective matrix ', num2str(i), ' out of ', num2str(ndb(2))]);
    
    % use each mode to create the matrix
    [N1_conv, N2_conv, ~] = c_rom.fe_model.get_OP_isaNavierStokes_conv(...
        c_rom.prim_shape(:,i), physics_opts);
    
    % project them onto the modes
    conv1_sys(:,:,i) = c_rom.prim_test'*N1_conv*c_rom.prim_shape;
    conv2_sys(:,:,i) = c_rom.prim_test'*N2_conv*c_rom.prim_shape;
end

% loop over the turbulence fields and get the diffusive parts
for i=1:ndf(2)
    % show the progress
    disp(['building diffusive matrix ', num2str(i), ' out of ', num2str(ndf(2))]);

    % change the turbulence field of the current stage
    physics_opts.mtvec = c_rom.visc_shape(:, i);
    [N_diff ,~] = c_rom.fe_model.get_OP_isaNavierStokes_diff(...
        zeros(ndb(1),1), physics_opts);

    % project them onto the modes
    diff_sys(:,:,i) = c_rom.prim_test'*N_diff*c_rom.prim_shape;
end
end

