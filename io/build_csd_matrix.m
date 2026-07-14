function [data_cov, chanlabels_out] = build_csd_matrix(freq_csd, chanlabels_wanted)
% Builds a chan x chan Hermitian CSD matrix from FieldTrip freq output,
% with explicit diagnostics for missing/NaN entries.
%
% freq_csd          : output of ft_freqanalysis (cfg.output = 'powandcsd')
% chanlabels_wanted : cell array, the canonical channel order to enforce
%                      (e.g. matching dics_result.avg.filter ordering)
% freqIdx           : index into crsspctrm's frequency dimension (1 if single foi)

Nchan = numel(chanlabels_wanted);
data_cov = nan(Nchan, Nchan);

% Map wanted labels to indices for fast lookup
[tf, locA] = ismember(chanlabels_wanted, chanlabels_wanted); %#ok<ASGLU>
labelIndex = containers.Map(chanlabels_wanted, 1:Nchan);

%% --- Step 1: fill off-diagonal / whatever pairs labelcmb gives us ---
nFilled = 0;
nSkippedNotWanted = 0;
nSkippedNoMatch = 0;

for k = 1:size(freq_csd.labelcmb, 1)
    lab1 = freq_csd.labelcmb{k,1};
    lab2 = freq_csd.labelcmb{k,2};

    if ~isKey(labelIndex, lab1) || ~isKey(labelIndex, lab2)
        nSkippedNotWanted = nSkippedNotWanted + 1;
        continue
    end

    i = labelIndex(lab1);
    j = labelIndex(lab2);

    val = freq_csd.crsspctrm(k, 1);

    data_cov(i,j) = val;
    data_cov(j,i) = conj(val);   % Hermitian mirror — safe even if i==j (val should be real then)

    nFilled = nFilled + 1;
end

fprintf('Filled %d entries from labelcmb (skipped %d not-in-chanlabels_wanted)\n', ...
    nFilled, nSkippedNotWanted);

%% --- Step 2: check whether diagonal got filled; if not, try powspctrm ---
diagIdx = 1:Nchan;
diagVals = data_cov(sub2ind(size(data_cov), diagIdx, diagIdx));
nNanDiag = sum(isnan(diagVals));

if nNanDiag > 0
    fprintf('%d/%d diagonal entries still NaN after labelcmb pass.\n', nNanDiag, Nchan);

    if isfield(freq_csd, 'powspctrm')
        fprintf('Attempting to fill diagonal from freq_csd.powspctrm...\n');
        for ci = 1:Nchan
            lab = chanlabels_wanted{ci};
            pIdx = find(strcmp(freq_csd.label, lab));
            if ~isempty(pIdx) && isnan(data_cov(ci,ci))
                data_cov(ci,ci) = freq_csd.powspctrm(pIdx, 1);
            end
        end
    else
        warning('No powspctrm field found — diagonal cannot be auto-recovered. Check cfg.output.');
    end
end

%% --- Step 3: final diagnostics ---
diagVals = data_cov(sub2ind(size(data_cov), diagIdx, diagIdx));
fprintf('Diagonal after recovery attempt: %d/%d still NaN, %d with nonzero imag part > 1e-10\n', ...
    sum(isnan(diagVals)), Nchan, sum(abs(imag(diagVals)) > 1e-10));

nNanTotal = sum(isnan(data_cov(:)));
if nNanTotal > 0
    fprintf('WARNING: %d/%d total matrix entries are still NaN.\n', nNanTotal, Nchan^2);
    [rIdx, cIdx] = find(isnan(data_cov));
    badPairs = unique([rIdx, cIdx], 'rows');
    fprintf('Example missing entries (first 5):\n');
    for r = 1:min(5, size(badPairs,1))
        fprintf('  (%s, %s)\n', chanlabels_wanted{badPairs(r,1)}, chanlabels_wanted{badPairs(r,2)});
    end
end

% Force exact Hermitian symmetry (kills float-level asymmetry)
data_cov = (data_cov + data_cov') / 2;

% Defensive: diagonal should be real
if any(abs(imag(diag(data_cov))) > 1e-10)
    warning('Diagonal has nontrivial imaginary part — check upstream computation.');
end
data_cov(sub2ind(size(data_cov),1:Nchan,1:Nchan)) = real(diag(data_cov));

chanlabels_out = chanlabels_wanted;
end