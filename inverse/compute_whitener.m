function [Cy, W] = compute_whitener(Cy)

[Un, Sn2] = svd(Cy, 'econ');
Sn = sqrt(diag(Sn2));
tol = length(Sn) * eps(single(Sn(1)));
Rank_Noise = sum(Sn > tol);

% Rebuild cov
Un = Un(:, 1:Rank_Noise);
Sn = Sn(1:Rank_Noise);
Cy = Un * diag(Sn.^2)*Un';

% Regularize (In this case by median Singluar value; like in LCMV from
% bstrm)
Sn = max(Sn, median(Sn));
% Rebuild cov again
Cy = Un * diag(Sn.^2) * Un';

% Build Whitener (In brainstorm the inverse whitener is build
%iW = Un * diag(1./Sn)*Un';
% I = W * iW
W = Un * diag(Sn) * Un';