function X = normball_samp(p,n,R,nvec)
%NORMBALL_SAMP generates samples within the p-Norm limit R
%   generates nvec samples of a n-dimensional parameter space within the
%   p-norm thershold defined by R.
%   if p is 2 then the norm is the euclidean one and the sampling will be
%   within a n-dimensional ball (a circle if n = 2) with radius R.
%
%   source: https://de.mathworks.com/matlabcentral/answers/530358-generate-a-vector-under-norm-condition

% if nvec was not provided, just generate one vector.
if nargin < 4
    nvec = 1;
end

% generate a uniform sample, of size nvec, that lies beRtween 0 and chi2cdf(n,r^2)
u = rand(1,nvec)*chi2cdf(R^2,n);
% invert those numbers through the chi-square CDF
% then take the sqrt. This uses the classic inverse CDF method to generate numbers
% from the desired distribution. Take the sqrt, since we really wanted a chi sample.
r = sqrt(chi2inv(u,n));

% generate an unconstrained sample in n dimensions via random sampling
X = randn(n,nvec); %X = lhsdesign(n, nvec);

% scale the samples so they ALL have exactly unit norm, then rescale them
% so the norm will be r.
xnorms = vecnorm(X,p,1);

% in case someone has an old MATLAB release than R2016b, use bsxfun
% X = bsxfun(@times,X,r./xnorms);
X = (X.*(r./xnorms))';

end

