function set_sys_mat(self, opts)
%SET_SYS_MAT Set the systems matrices for the reduced model
%   Using the current base functions, the systems matrices are calculated
%   and stored in root and as properties of the model

% construct the path to the folder and sys file
sys_folder = fullfile(self.root_folder, self.base_name);
sys_file = fullfile(sys_folder, 'sys_matrices.mat');

% get the 3rd order tensors for the reduced model
if isfile(sys_file)
    disp('loading stored system matrices...')
    sys_matrices = load(sys_file).sys_matrices;
elseif ~isfile(sys_file)
    [diff_sys, conv_sys] = self.get_sys_build(opts);
    sys_matrices.diff_sys = diff_sys;
    sys_matrices.conv_sys = conv_sys;
    save(sys_file, 'sys_matrices');
end

% set the matrices as properties for the rom model
self.diff_sys = sys_matrices.diff_sys;
self.conv_sys = sys_matrices.conv_sys;
end

