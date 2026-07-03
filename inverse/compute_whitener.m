function iW = compute_whitener(Cy, params)

method = 'none';

if exist('params', 'var')
    if isfield(params, 'method')
        method = params.method;
    end
    if isfield(params, 'NoiseReg')
        NoiseReg = params.NoiseReg;
    end
end

[Un, Sn2] = svd(Cy, 'econ');
Sn = sqrt(diag(Sn2));
tol = length(Sn) * eps(single(Sn(1)));
Rank_Noise = sum(Sn > tol);

% Rebuild cov
Un = Un(:, 1:Rank_Noise);
Sn = Sn(1:Rank_Noise);


% Bild whitener. s^-0.5 = 1/sqrt(s) -> before we took already the sqrt(s),
% so now computing the whitener as u * (1./s) * u' = u * s^-0.5 * u'!!
if strcmp(method, 'none')
    iW = Un * diag(1./Sn) * Un';
elseif strcmp(method, 'median')
    % Regularize (In this case by median Singluar value; like in LCMV from
    % bstrm)
    Sn = max(Sn, median(Sn));
    
    % Build inverse whitener 
    iW = Un * diag(1./Sn) * Un';
    % W = Un * diag(Sn) * Un';
    % I = W * iW
elseif strcmp(method, 'reg')
    if ~exist('NoiseReg', 'var')
        error('For the method "reg", please specify a NoiseReg parameter!');
    end

    RidgeFactor = mean(diag(Sn2)) * NoiseReg; % Hamalainen's preferred measure

    iW = Un*diag(1./sqrt(Sn.^2 + RidgeFactor))*Un'; % inverse whitener, symmetric
end

