function th = tanh_dist(c,n, varargin)
%TANH_DIST get a tanh distribution from 0 to 1

% create the distribution
th = (tanh(linspace(-c,c,n)) + 1)/2;
th = (th - th(1))/(th(end) - th(1));

% map it onto the interval (if given)
if ~isempty(varargin)
    mi = min(varargin{1}); ma = max(varargin{1});
    th = th*(ma - mi) + mi;
end
end

