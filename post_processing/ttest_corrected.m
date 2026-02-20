function [tval,p] = ttest_corrected(x, varargin)

% Set defaults
defaults = struct('correction', 0, 'pop_mean_diff', 0, 'tail', 0);
% Parse inputs
params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
    end
end

% Process remaining arguments
dim = find(size(x) ~= 1, 1);
if isempty(dim), dim = 1; end

nans = isnan(x);
if any(nans(:))
    samplesize = sum(~nans,dim);
else
    samplesize = size(x,dim); % a scalar, => a scalar call to tinv
end
df = max(samplesize - 1,0);
xmean = mean(x,dim, 'omitmissing');
varpop = var(x,[],dim, 'omitmissing');
if params.correction==0
    params.correction=.01*max(varpop);
end
corrsdpop=sqrt(varpop+params.correction);

ser = corrsdpop ./ sqrt(samplesize);
tval = (xmean - params.pop_mean_diff) ./ ser;
% Compute the correct p-value for the test, and confidence intervals
% if requested.
if params.tail == 0 % two-tailed test
    p = 2 * tcdf(-abs(tval), df);
elseif params.tail == 1 % right one-tailed test
    p = tcdf(-tval, df);
elseif params.tail == -1 % left one-tailed test
    p = tcdf(tval, df);
end
