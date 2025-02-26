function opts = get_opts(varargin)
%GET_OPTS setup the options for the solver
%   the solver needs various options like damping, initialization and so
%   on. this fuction provides the necessary struct, filling not used
%   options with default values.

% parse the inputs and get the opts
ip = inputParser();
ip.addParameter('default_case', 1);
ip.addParameter('primary_bound_th', 0.9);
ip.addParameter('primary_domain_th', 0.9);
ip.addParameter('viscos_bound_th',0.9);
ip.addParameter('viscos_domain_th',0.9);
ip.addParameter('geo_variations', {});
ip.addParameter('n_geo_diff', 4);
ip.addParameter('n_geo_conv', 4);
ip.addParameter('shape2test', 2);
ip.addParameter('build', 'approx')
ip.addParameter('R', 287);
ip.addParameter('kappa', 1.4);
ip.addParameter('ref_pressure', 101300);
ip.addParameter('ref_temperature', 298);
ip.addParameter('alpha', 1);
ip.addParameter('plot', true);
ip.parse(varargin{:});

% translate the inputs into a opts struct
opts = ip.Results;

end

