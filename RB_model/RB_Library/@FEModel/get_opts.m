function opts = get_opts(varargin)
%GET_OPTS setup the options for the solver
%   the solver needs various options like damping, initialization and so
%   on. this fuction provides the necessary struct, filling not used
%   options with default values.

% parse the inputs and get the opts
ip = inputParser();
ip.addParameter('damping', 1);
ip.addParameter('init', nan);
ip.addParameter('stab',1);
ip.addParameter('lin', 1);
ip.addParameter('mue', 1.72e-5);
ip.addParameter('rho', 1.225);
ip.addParameter('plot', true);
ip.addParameter('time_rec', false);
ip.addParameter('gamma', 1.4);
ip.parse(varargin{:});

% translate the inputs into a opts struct
opts = ip.Results;

end

