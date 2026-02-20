function compare_power_change_subject(sFiles, varargin)

% Set defaults
defaults = struct('roi_cutoff', 60, 'tail', 0, 'use_scouts', 0);
% Parse inputs
params = struct()
disp('Testing still active. If you want to run this script from another script disable testing!');
%params = struct(varargin{:});
for f = fieldnames(defaults)',
    if ~isfield(params, f{1}),
        params.(f{1}) = defaults.(f{1});
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

    % Compute global roi
    %% To implement: make thresholding param variable based on user input, make other options of defining ROIs possible e.g. via Scouts!
    if ~params.use_scouts
        pial_avg = mean(pial_diff(:, 1, i), 2);
        pial_thresh = prctile(pial_avg, params.roi_cutoff);
        pial_mask = pial_avg > pial_thresh;
        % Use Leadfield
        white_avg = mean(white_diff(:, 1, i), 2);
        white_thresh = prctile(white_avg, params.roi_cutoff);
        white_mask = white_avg > white_thresh;
    else
        % Do something
    end
    
    
    multilayer_mask = pial_mask | white_mask;
    fprintf('%i vertices (%.2f%s of all vertices) are included in comparison.\n', sum(multilayer_mask), sum(multilayer_mask)/length(multilayer_mask), '%');

    % Average over over vertices in ROI
    avg_trial_change = mean(pial_white_diff(multilayer_mask, :, i));
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
end

disp(t_vals);
disp(p_vals);

