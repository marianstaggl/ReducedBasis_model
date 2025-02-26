function base_fun = get_base_isa_naive(self, opts)
%GET_BASE_FUN read data from root and create the base functions
%   The solution data is read from the root folder and the fields are
%   decomposed into multiple modes.

% read the field of primary variables and the turbulence field
all_primary = load(fullfile(self.root_folder, 'all_primary')).all_primary;
if ~isequal(all_primary.variables, self.fe_model.sys_variables)
    error(['the model only works with systems variables ',...
        num2str(self.fe_model.sys_variables)]);
end
all_viscos = load(fullfile(self.root_folder, 'all_viscos')).all_viscos;

% identify the default solution and center the data around
mean_prim = all_primary.values(:, opts.default_case);
prim_cent = all_primary.values - mean_prim;
mean_visc = all_viscos.values(:, opts.default_case);
visc_cent = all_viscos.values - mean_visc;

% decompose the primary and viscos fields
[bound_prim, bp_sel, bp_decay] = decompose_field(prim_cent, ...
    opts.primary_bound_th, self.bound_nodes.primary);
[domain_prim, dp_sel, dp_decay] = decompose_field(bound_prim(:, ~bp_sel), ...
    opts.primary_domain_th);
[bound_visc, bv_sel, bv_decay] = decompose_field(visc_cent, ...
    opts.viscos_bound_th); % , self.bound_nodes.viscos

% build the full shape and test function base and the one for viscos
n_shape = sum(dp_sel); 
base_fun.prim_shape = [mean_prim, bound_prim(:, bp_sel), ...
    domain_prim(:, 1:n_shape)];
n_test = min([4*n_shape, size(domain_prim, 2)]); % allow a max of 4 for overdetermination
base_fun.prim_test = [mean_prim, bound_prim(:, bp_sel), ...
    domain_prim(:, 1:n_test)];
base_fun.visc_shape = [mean_visc, bound_visc(:, bv_sel)];

% get the selector for the boundary modes
base_fun.sys_eqn = 'isa. navier stokes';
base_fun.sys_variables = all_primary.variables;
base_fun.prim_bound = false(1, size(base_fun.prim_shape, 2));
base_fun.prim_bound(1:sum(bp_sel) + 1) = true;
base_fun.visc_bound = false(1, size(base_fun.visc_shape, 2));
base_fun.visc_bound(1:sum(bv_sel) + 1) = true;

if opts.plot, plot_decay(bp_decay, dp_decay, bv_decay); end
end

function plot_decay(prim_bc, prim_do, visc_bc)
subplot(2, 1, 1);
semilogy(prim_bc); hold on; semilogy(prim_do);
legend({'boundary modes', 'domain modes'}); grid on;
xlabel('# mode'); ylabel('explained energy')
subplot(2, 1, 2);
semilogy(visc_bc); hold on; grid on;
legend({'viscosity field'}); xlabel('# mode'); 
ylabel('explained variance')
end

function [modes, mode_sel, singular_values] = decompose_field(field, th,...
    varargin)
% check if all nodes are subject to decomposition
if isempty(varargin)
    [modes, s, ~] = svd(field, 'econ');
else
    [~, s, v] = svd(field(varargin{1}, :));
    modes = field * v;
end

% calculate the energy containing level
singular_values = diag(s);
if th == 0 % in this case the model is a gappyPOD model
    mode_sel = false(1, size(field,2));
elseif th > 0 && th < 1
    e_level = cumsum(diag(s)) / sum(s, 'all');
    e_level = [0, reshape(e_level(1:end-1), 1, [])];
    mode_sel = ~(e_level > th);
elseif th >= 1
    mode_sel = false(1, size(field,2));
    mode_sel(1:th) = true;
else
    error("this shouldn't be possible")
end
end

