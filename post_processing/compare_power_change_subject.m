function compare_power_change_subject(sFiles, params)

% Set defaults
roi_cutoff = 60;
tail = 0;
roi = [];
verbose = true;

% Parse inputs
if exist('params', 'var') && ~isempty(params)
    if isfield(params, 'roi_cutoff')
        roi_cutoff = params.roi_cutoff;
    end

    if isfield(params, 'tail')
        tail = params.tail;
    end

    if isfield(params, 'roi')
        roi = params.roi;
    end

    if isfield(params, 'verbose')
        verbose = params.verbose;
    end
end

sData = in_bst_data(sFiles{1});
n_contrasts = size(sData.Freqs, 1);

surf = in_bst_data(sData.SurfaceFile);
n_Vertices = size(surf.Vertices, 1);
vertices_white = 1:n_Vertices/2;
vertices_pial = n_Vertices/2+1:n_Vertices;

pial_diff = zeros(size(surf.Vertices, 1)/2, length(sFiles), n_contrasts);
white_diff = zeros(size(surf.Vertices, 1)/2, length(sFiles), n_contrasts);
pial_white_diff = zeros(size(surf.Vertices, 1)/2, length(sFiles), n_contrasts);

%% Load all trial differences for all contrasts
for i = 1:length(sFiles)
    sTrial = in_bst_data(sFiles{i});
    for ii=1:n_contrasts
        % Split into white and pial surface results
        white_diff(:, 1, ii) = sTrial.TF(vertices_white, 1, ii);
        pial_diff(:, 1, ii) = sTrial.TF(vertices_pial, 1, ii);
        % calculate difference between the two surfaces
        % Here taking the absolute is very! important
        pial_white_diff(:, i, ii) = abs(pial_diff(:, 1, ii)) - abs(white_diff(:, 1, ii));
    end
end

contrast_name = {};
t_vals = zeros(n_contrasts, 1);
p_vals = zeros(n_contrasts, 1);

for i=1:n_contrasts

    if isempty(roi)
        % Compute global roi
        pial_avg = mean(pial_diff(:, 1, i), 2);
        pial_thresh = prctile(pial_avg, roi_cutoff);
        pial_mask = pial_avg > pial_thresh;
        % Use Leadfield
        white_avg = mean(white_diff(:, 1, i), 2);
        white_thresh = prctile(white_avg, roi_cutoff);
        white_mask = white_avg > white_thresh;
    else
        % Compute roi based on specified scout
        pial_mask = zeros(size(pial_diff, 1), 1);
        pial_mask(roi) = true; 
        % Use Leadfield
        white_mask = zeros(size(pial_diff, 1), 1);
        white_mask(roi) = true;    
    end
    
    multilayer_mask = pial_mask | white_mask;
    fprintf('%i vertices (%.2f%s of all vertices) are included in comparison.\n', sum(multilayer_mask), 100*(sum(multilayer_mask)/length(multilayer_mask)), '%');

    % Average over over vertices in ROI
    if sum(multilayer_mask) == 1
        avg_trial_change = pial_white_diff(multilayer_mask, :, i);
    else
        avg_trial_change = mean(pial_white_diff(multilayer_mask, :, i), 2);
    end    
    
    contrast_name{i} = sTrial.Freqs{i};

    % compare vertex activity against 0 -> adapted std of population
    % (https://doi.org/10.1016/j.neuroimage.2011.10.027)
    % Here in bonaiutos code the correction param was choosen to be
    % 25*var(avg_trial_change) -> does not go with the explanation in above
    % paper does it? Default var scaling in ttest_corrected is by 0.01
    % (also 10*higher then recommended by above paper
    % Nachdenken!
    [t_vals(i), p_vals(i)] = ttest_corrected(avg_trial_change, 'correction', 10e-4*max(var(pial_white_diff(multilayer_mask, :, ii), [], 2)), ...
        'tail', 0);

    if verbose
        fprintf('Condition: %s\n', contrast_name{i});
        fprintf('t-value: %.3f, p-value: %.3f\n\n', t_vals(i), p_vals(i));
    end
end