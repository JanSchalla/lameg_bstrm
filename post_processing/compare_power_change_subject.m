function [contrast_name, t_vals, p_vals] = compare_power_change_subject(sFiles, params)

if ~brainstorm('status')
    brainstorm nogui
end

% Set defaults
roi_cutoff = 70;
tail = 0;
roi = [];
verbose = true;
create_vis = true;

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

    if isfield(params, 'create_vis')
        create_vis = params.create_vis;
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
        white_diff(:, i, ii) = sTrial.TF(vertices_white, 1, ii);
        pial_diff(:, i, ii) = sTrial.TF(vertices_pial, 1, ii);
        % calculate difference between the two surfaces
        % Here taking the absolute is very important!
        pial_white_diff(:, i, ii) = abs(pial_diff(:, i, ii)) - abs(white_diff(:, i, ii));
    end
end

contrast_name = {};
t_vals = zeros(n_contrasts, 1);
p_vals = zeros(n_contrasts, 1);

for i=1:n_contrasts

    % Calcualte correction metric (https://doi.org/10.1016/j.neuroimage.2011.10.027)
    % Correction is applied below
    var_full_map = var(pial_white_diff(:, :, i), [], 2);
    delta = 1e-3 * max(var_full_map);

    if isempty(roi)
        % Compute global roi
        pial_t_statistic = ttest_corrected(pial_diff(:, :, i)');
        white_t_statistic = ttest_corrected(white_diff(:, :, i)');
        pial_thresh = prctile(pial_t_statistic, roi_cutoff);
        pial_mask = pial_t_statistic > pial_thresh;
        white_thresh = prctile(white_t_statistic, roi_cutoff);
        white_mask = white_t_statistic > white_thresh;
    else
        % Compute roi based on specified scout
        pial_mask = zeros(size(pial_diff, 1), 1);
        pial_mask(roi) = true; 
        white_mask = zeros(size(pial_diff, 1), 1);
        white_mask(roi) = true;    

        pial_t_statistic = ttest_corrected(pial_diff(logical(pial_mask), :, i)');
        white_t_statistic = ttest_corrected(white_diff(logical(white_mask), :, i)');
    end

    if create_vis

        if length(split(sFiles{1}, "/")) > length(split(sFiles{1}, "\"))
            separator = "/";
        else
            separator = "\";
        end
        
        pial_white_diff_t_statistic = ttest_corrected(pial_white_diff(:, :, i)');
        pial_white_diff_var         = var (pial_white_diff(:, :, i)');
        pial_white_diff_mean        = mean(pial_white_diff(:, :, i)');
        
        tokens = split(sFiles{1}, separator);
        sSubject = bst_get('Subject', tokens{1});
        single_layer_surf_idx = find(contains({sSubject.Surface.Comment}, 'corresponding'));
        
        nSrc = round(n_Vertices / 2);
        
        tf_template = db_template('timefreq');
        
        tf_template.Comment     = sprintf('Multilayer Comparison Results | Condition: %s', sTrial.Freqs{i});
        tf_template.SurfaceFile = sSubject.Surface(single_layer_surf_idx(1)).FileName;
        tf_template.HeadModelFile = ''; 
        tf_template.DataFile      = '';
        tf_template.DataType      = 'results';
        tf_template.RowNames      = cellfun(@num2str, num2cell(1:nSrc), 'UniformOutput', false);
        tf_template.Time = [0, 0];   
        tf_template.Measure       = 'power';
        tf_template.Method        = 'hilbert';
        tf_template.nAvg          = 1;
        tf_template.Leff          = 1;
        
        % TF is [nSources × nTime × nFreqs]
        tf_template.TF = zeros(nSrc, 1, 5);
        tf_template.TF(:, 1, 1) = pial_t_statistic;
        tf_template.TF(:, 1, 2) = white_t_statistic;
        tf_template.TF(:, 1, 3) = pial_white_diff_t_statistic;
        tf_template.TF(:, 1, 4) = pial_white_diff_mean;
        tf_template.TF(:, 1, 5) = pial_white_diff_var;
        
        % Freqs cell: {label, fmin, fmax} — the label column is what appears in the GUI slider
        tf_template.Freqs = {
            't (Pial)',             sTrial.Freqs{i,2}, sTrial.Freqs{i,3};
            't (White)',            sTrial.Freqs{i,2}, sTrial.Freqs{i,3};
            't (P-W)',  sTrial.Freqs{i,2}, sTrial.Freqs{i,3};
            'Diff Mean (P-W)',         sTrial.Freqs{i,2}, sTrial.Freqs{i,3};
            'Diff Var (P-W)',     sTrial.Freqs{i,2}, sTrial.Freqs{i,3};
        };
        
        tf_template.History = {
            datestr(now,'dd/mm/yy-HH:MM'), 'compute', 'Save multilayer whole brain results'
        };
        
        % BugFix Visualization (Dirty Fixc, but working for the moment)
        tf_template.GridLoc    = [];
        tf_template.GridAtlas  = [];

        intra_path = fileparts(file_fullpath(sSubject.FileName));
        intra_path = replace(intra_path, 'anat', 'data');
        
        out_fname = sprintf('timefreq_psd_multilayer_result_%s.mat', sTrial.Freqs{i});
        save(fullfile(intra_path, '@intra', out_fname), '-struct', 'tf_template');

    end

    multilayer_mask = pial_mask | white_mask;
    
    if verbose
        fprintf('%i vertices (%.2f%s of all vertices) are included in comparison.\n', sum(multilayer_mask), 100*(sum(multilayer_mask)/length(multilayer_mask)), '%');
    end

    % Average over vertices in ROI

    avg_trial_change = mean(pial_white_diff(multilayer_mask, :, i), 1);  
    
    contrast_name{i} = sTrial.Freqs{i};
    
    [t_vals(i), p_vals(i)] = ttest_corrected(avg_trial_change, 'correction', delta, ...
        'tail', 0);

    if verbose
        fprintf('Condition: %s\n', contrast_name{i});
        fprintf('t-value: %.3f, p-value: %.3f\n\n', t_vals(i), p_vals(i));
    end
end