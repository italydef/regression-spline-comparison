function tdata = randt(dim, nvecs, mean, varcov)
%
% randt(dim, nvecs, mean, varcov)
%    returns a matrix of dimension dim by nvecs of observations
%    from a multivariate (dim dimension) t distribution with the provided
%    mean and varcov on 4 degrees of freedom.
%    The columns are the multivariate t r.v.s
%
%    Written by Elizabeth Slate
%

if nargin < 4
  varcov = eye(dim);
end
if nargin < 3
  mean = zeros(dim,1);
end

df = 4;

% generate the appropriate number of standard Gaussians
X = randn(dim, nvecs);

% generate nvecs gamma random variables with parameters df/2 and df/2
S2 = randgam(1, nvecs, df/2, df/2);

S = S2 .^ 0.5;

stndT = X ./ (ones(dim,1) * S);

tdata = sqrtm(varcov) * stndT + mean * ones(1,nvecs);

end



