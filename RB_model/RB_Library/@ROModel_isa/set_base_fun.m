function set_base_fun(self, base_name, varargin)
%SET_BASE_FUN Set the base functions for the reduced model
%   use existing cfd solutions to extract the modes for the reduced model.
%   the modes can be specified via varargin or they are calculated using
%   default options. If they already exist they are loaded from root
%   instead of creating them again

ip = inputParser();
ip.addParameter('base_fun', []);
ip.parse(varargin{:});

% if the base already exists, load it from root
base_folder = fullfile(self.root_folder, base_name);
base_file = fullfile(base_folder, 'base_fun.mat');

% get the base functions for the reduced model
if isfile(base_file)
    disp('loading stored base...')
    base_fun = load(base_file).base_fun;
elseif ~isfile(base_file) && ~isempty(ip.Results.base_fun)
    mkdir_safe(base_folder);
    base_fun = ip.Results.base_fun;
    save(base_file, 'base_fun');
else
    error("neither a base function struct nor a valid folder is given...")
end

% set the base function as property of the object
self.base_name = base_name;
self.prim_shape = base_fun.prim_shape;
self.prim_test = base_fun.prim_test;
self.prim_bound = base_fun.prim_bound;
self.visc_shape = base_fun.visc_shape;
self.visc_bound = base_fun.visc_bound;

% update the finite element model equations
self.fe_model.sys_equations = base_fun.sys_eqn;
self.fe_model.sys_variables = base_fun.sys_variables;

% delete the system matrices
self.diff_sys = [];
self.conv_sys = [];
end

function mkdir_safe(new_dir)
% only make the directory if it doesn't exist yet
if not(isfolder(new_dir))
    mkdir(new_dir)
end
end