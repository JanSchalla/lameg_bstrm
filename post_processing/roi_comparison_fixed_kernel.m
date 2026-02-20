function roi_comparison_fixed_kernel(sFiles)

template_session_results = load('/home/jan/Jan/993_TEMPLATES/matlab/bstrm_structures/template_multilayer_desc_results.mat');
template_inference_results = load('/home/jan/Jan/993_TEMPLATES/matlab/bstrm_structures/template_multilayer_stat_results.mat');

% sFiles cell of relaitves pathes to the trial data from the brainstorm db
sData = in_bst_data(sFiles{1})
time = sData.Time;
clear sData
%% Check if all files are of the same run
nTrialsTotal = length(sFiles);
ses_ids = unique(cellfun(@fileparts, sFiles, 'UniformOutput', false));

session_files = cell(1, length(ses_ids));
for ses = 1:length(ses_ids)
    file_mask = cellfun(@(x) contains(x, ses_ids{ses}), sFiles);
    session_files{ses} = {sFiles{file_mask}};
end

fprintf('%i Sessions identfied.\n', length(ses_ids));

%% 1st: Let user define baseline, window of interest (woi) and frequency of interest (foi).
prompt = {'Specify baseline time window (s, e. g. -0.5, -0.1):', ...
    'Specify window of interest (s, e. g. 0, 1):', ...
    'Specify frequency of interest (Hz, e. g. 13, 30)'};
dlgtitle = 'Surface Power Comparison';
fieldsize = [1 45; 1 45; 1 45];
definput = {'','', ''};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

% Parse user inputs
win_base = str2double(split(answer{1}, ','));
win_oi = str2double(split(answer{2}, ','));
freq_oi = str2double(split(answer{3}, ','));

% Check if trial bounds are valid
if min(win_base) < min(time)
    error('Baseline starts outside of trial. Baseline cannot be defined smaller then %.3f seconds', min(time));
end

if max(win_base) > min(win_oi)
    warning('Baseline and Window if Interest overlap. This is not recommended.');
    dlgTitle    = 'Overlap of Baseline and WOI';
    dlgQuestion = 'Do you wish to continue?';
    choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'No');

    if strcmp(choice, 'No')
        error('Baseline and WOI overlap.');
    end
end

if max(win_oi) > max(time)
    error('Window of Interest ends outside of trial bound. WOI cannot be defined bigger then %.3f seconds', max(time));
end

% Check if frequency resolution is enough
sfreq = round(1/(time(2) - time(1)));

if max(freq_oi) > sfreq/2
    error('Selected Frequency of Interest is bigger then Nyquist Frequency (%i Hz).', round(sfreq/2));
end

% bring time information in sample space
[~, base_min_idx] = min(abs(time - win_base(1)));
[~, base_max_idx] = min(abs(time - win_base(2)));
win_base = [base_min_idx base_max_idx];
clear('base_min_idx', "base_max_idx");

[~, woi_min_idx] = min(abs(time - win_oi(1)));
[~, woi_max_idx] = min(abs(time - win_oi(2)));
win_oi = [woi_min_idx woi_max_idx];
clear('woi_min_idx', "woi_max_idx");
%% Load in Multilayer Kernel
%sKernel = in_bst_data('/home/jan/Jan/brainstorm_db/00_layerAnalysis_Test/data/pilot001/RDKPILOT001_run03_20250117_(3)_notch_resample/results_PNAI_MEG_GRAD_MEG_MAG_KERNEL_251215_1443.mat');
used_SR = [];
pial_baseline_all = [];
white_baseline_all = [];
pial_vals_all = [];
white_vals_all = [];

for ses = 1:length(ses_ids)
    % init brainstorm
    bstrm_path = session_files{ses}{1};
    if length(split(bstrm_path, "/")) > length(split(bstrm_path, "\"))
        separator = "/";
    else
        separator = "\";
    end
    tokens = split(bstrm_path, separator);
    protocol_id = find(ismember(tokens, "brainstorm_db"), 1) + 1;
    % chekc if previous path was relative
    if isempty(protocol_id)
        bstrm_path = file_fullpath(bstrm_path);
        % Redo if previous path was relative
        tokens = split(bstrm_path, separator);
        protocol_id = find(ismember(tokens, "brainstorm_db"), 1) + 1;
    end
    protocol_name = tokens{protocol_id};
    subject_id = find(ismember(tokens, "data"), 1) + 1;
    subject_name = convertStringsToChars(tokens{subject_id});
    study_name = convertStringsToChars(tokens{subject_id + 1});
    data_name = convertStringsToChars(tokens{subject_id + 2});
    
    bstrm_out_path = fileparts(bstrm_path);
    % find sr file
    kernel_file = dir(fullfile(strjoin(tokens(1:subject_id+1), separator), '*KERNEL*.mat'));
    if isempty(kernel_file)
        fprintf('No source reconstruction found for %s. Cancelling process ...', subject_name);
        return
    else
        if length(kernel_file) == 1
            kernel_file = string(fullfile(kernel_file.folder, kernel_file.name));
            used_SR(ses) = kernel_file;
        else
            % Prepare list of full file paths and display names
            file_paths = fullfile({kernel_file.folder}, {kernel_file.name});
            file_names = {kernel_file.name};
            
            [selection, is_ok] = listdlg('PromptString', sprintf('Multiple KERNEL files found for %s.\nSelect one:', subject_name), ...
                                         'SelectionMode', 'single', ...
                                         'ListString', file_names, ...
                                         'Name', 'Choose Kernel File');
            
            if is_ok && ~isempty(selection)
                kernel_file = string(file_paths{selection});
                used_SR(ses) = kernel_file;
                fprintf('Selected kernel file: %s\n', kernel_file);
            else
                fprintf('No kernel file selected for %s. Cancelling process ...\n', subject_name);
                return
            end
        end
    end

    sKernel = load(kernel_file);
    
    nSources = size(sKernel.ImagingKernel, 1);
    nVertices = nSources/2;
    
    nSources_white = 1:nVertices;
    nSources_pial = nVertices+1:nSources;
    
    if length(nSources_pial) ~= length(nSources_white)
        error('White and pial subsurfaces contain different nr. of vertices! Aborting...');
    end

    sesTrials = length(session_files{ses});
    
    pial_baseline_ses = zeros(nVertices, sesTrials);
    white_baseline_ses = zeros(nVertices, sesTrials);
    pial_vals_ses = zeros(nVertices, sesTrials);
    white_vals_ses = zeros(nVertices, sesTrials);
    f = waitbar(0, 'Extracting surface specific power');
    for iTrial = 1:sesTrials
        waitbar(iTrial/sesTrials, f, sprintf('Extracting surface specific power: %d %%', floor(iTrial/sesTrials*100)));
        sRawData = in_bst_data(session_files{ses}{iTrial});
    
        sensorData = sRawData.F(sKernel.GoodChannel, :);
        
        % Here one can adapt in which order project the sources and filter.
        % The kernel however is computed on the unfiltered sensor data
        % Commented out for the sake of local testing
        % [ebase_amp, ewoi_amp] = SrFilterPowerExtract(sensorData, sKernel, win_base, win_of, freq_of, sfreq);
        
        %% Testing block (remove when test works)
        %% Testing SrFilterPowerExtract
        % First project the data into source space, then filter each source
        source_data = sKernel.ImagingKernel * sensorData;
        source_data_bp = bandpass(source_data', [freq_oi(1) freq_oi(2)], sfreq);

        %% Here i can take the hilbert envelope and take the mean in the specified timewindow
        % Or should i use pwelch or timefrquency decomposition again? Same as
        % above

        source_data_power = abs(hilbert(source_data_bp)).^2;

        base_power = mean(source_data_power(win_base(1):win_base(2), :)', 2);
        woi_power = mean(source_data_power(win_oi(1):win_oi(2), :)', 2);
    
        %% Testing SrPowerExtract 
        % source_data = sKernel.ImagingKernel * sensorData;
        % 
        % % demean each vertex
        % source_data_demean = source_data - mean(source_data, 2);
        % 
        % % Compute TF repreentation [Channels/Vertices x timepoints x freqs]
        % P = morlet_transform(source_data_demean, time, freq_of(1):freq_of(2), 1, 3, 'y');
        % 
        % % Extract mean power for the baseline and woi for each vertex
        % base_power = mean(P(: ,win_base(1):win_base(2), :), 2:3);
        % woi_power = mean(P(:, win_of(1):win_of(2), :), 2:3);
        % 
        % % clean up
        % clear P
        %%
    
        white_baseline_ses(:, iTrial) = base_power(nSources_white);
        pial_baseline_ses(:, iTrial) = base_power(nSources_pial);
    
        white_vals_ses(:, iTrial) = woi_power(nSources_white);
        pial_vals_ses(:, iTrial) = woi_power(nSources_pial);

        % if SaveIndividualTrial
        %     % save individual trials
        % end
    end
    close(f);

    %test = in_bst_data(sKernel.HeadModelFile); % the forward field is in 3D (x, y, z), how can i only take the top 60%? Combine the forward field? 
    ses_template = template_session_results;
    ses_template.TF = zeros(nVertices, 1, 2);
    ses_template.HeadModelFile = sKernel.HeadModelFile;
    ses_template.SurfaceFile = sKernel.SurfaceFile; % This i need to adapt or i just take the results double??
    ses_template.Comment = ses_ids(ses);
    ses_template.Time = [min(time) max(time)];
    ses_template.Method = 'hilbert';
    % Here we calculate the change of activity from the baseline
    %% Should this be saved?
    white_diff_ses = white_vals_ses - white_baseline_ses;
    pial_diff_ses= pial_vals_ses - pial_baseline_ses;
    % Do we need the stats for the separate surfaces?
    %[H,pvals,ci,STATS]=ttest(pial_diff');
    
    % Here we substract the absolute change of activity on the pial surface
    % from the white matter surface
    pial_white_diff_ses= abs(pial_diff_ses) - abs(white_diff_ses);
    
    % With ttest we test if the distribution is different from 0
    [H, pvals, ci, STATS] = ttest(pial_white_diff_ses');
    % Here we extract the t-values for each vertex
    pial_white_tvals_ses=STATS.tstat';

    ses_template.TF(:, 1, 1)= pial_white_tvals_ses;
    ses_template.Freqs{1, 1} = 'beta';
    ses_template.Freqs{1, 2} = sprintf('%i, %i', freq_oi(1), freq_oi(2));
    ses_template.Freqs{1, 3} = 'Pial - White t-values';

    ses_template.TF(:, 1, 2)= pvals;
    ses_template.Freqs{2, 1} = 'beta';
    ses_template.Freqs{2, 2} = sprintf('%i, %i', freq_oi(1), freq_oi(2));
    ses_template.Freqs{2, 3} = 'Pial - White p-values';

    out_fname = 'timefreq_psd_average_multilayer.mat';
    save(fullfile(bstrm_out_path, out_fname), '-struct', 'ses_template');

    % Here i average over trials to get the mean change per vertex
    average_pial_white_diff = mean(pial_white_diff_ses, 2);

    % Inside i average over verticies to get the mean change over trials in
    % the specified region!
    [t, p] = ttest_corrected(average_pial_white_diff, 'correction', 25*var(average_pial_white_diff))
    

    %% To do
    % Save data for each session

end

%% To do
% ----
% Combine session and save in intra subject folder
% ----
% ----
% write t-test function to comparte t-values/effect size/difference against 0
% To test between surfaces:
% dof for each session is nTrials-1
% ----
% The t-value is calculated by substracting the value from the mean and
% divide it by the std_error, also a p_val is calculated
% ----
% The data from the mask is then put into the ttest_corrected function
% there we estimate the population mean and population variance which is
% corrected by bonaiuto by adding 25*var(x) to get the corrected std -> why 25??
% with the corrected std the std_error is calcualted by dividing it with
% the sqrt(SampleSize)
% ----
% ----
% This can be done in the ttest
% I can based on the forward field apply a threshold on which vertices i
% take or not take (e.g. only take the 60% of vertices with the highest
% forward field activity (to test!)
% ----
% ----
% Decide with what power estimation method to extract the power values ->
% and/or test -> Test how it works with timefreq (very computation heavy),
% hampel and fft/pwelch
% ----
% ----
% Apply masking option when running the ttest
% Depending if i want to test it hemisphere specific, whole brain, ROI i
% need a mask -> this could be done using the bstrm atlases
% ----
% ----
%% Helper Functions
function [base_power, woi_power] = SrPowerExtract(sensor_data, time_vec, imgKernel, baseline, woi, foi, sfreq)

    %% Computationally very heavy (about 13.5s for the beta band (13-30Hz)
    % First project the data into the source space, and then extract the
    % frequency of interest from the defined time window.
    source_data = imgKernel * sensor_data;
    
    % demean

    % Compute TF repreentation [Channels/Vertices x timepoints x freqs]
    P = morlet_transform(source_data, time_vec, foi(1):foi(2), 1, 3, 'y');

    % Extract mean power for the baseline and woi for each vertex
    base_power = mean(P(: ,baseline(1):baseline(2), :), 2:3);
    woi_power = mean(P(:, woi(1):woi(2), :), 2:3);
end

function [base_amp, woi_amp] = SrFilterPowerExtract(sensor_data, imgKernel, baseline, woi, foi, sfreq)
    % First project the data into source space, then filter each source
    source_data = imgKernel * sensor_data;
    source_data_bp = bandpass(source_data', [foi(1) foi(2)], sfreq);

    %% Here i can take the "instantaneous" power and take the mean in the specified timewindow
    % Or should i use pwelch or timefrquency decomposition again? Same as
    % above

    source_data_power = abs(hilbert(source_data_bp)).^2;
    
    base_amp = mean(source_data_power(baseline(1):baseline(2), :)', 2);
    woi_amp = mean(source_data_power(woi(1):woi(2), :)', 2);
end

function [base_amp, woi_amp] = FilterSrPowerExtract(sensor_data, imgKernel, baseline, woi, foi, sfreq)
    % First filter the data in the sensor space, then project them into
    % source space
    sensor_data_bp = bandpass(sensor_data', [foi(1) foi(2)], sfreq);
    source_data = imgKernel * sensor_data_bp';

    %% Here i can take the hilbert envelope and take the mean in the specified timewindow
    % Or should i use pwelch or timefrquency decomposition again? Same as
    % above
    source_data_power = abs(hilbert(source_data_bp)).^2;
    
    base_amp = mean(source_data_power(baseline(1):baseline(2), :)', 2);
    woi_amp = mean(source_data_power(woi(1):woi(2), :)', 2);
end

% Example usages
%pwelch(LFP_right_bp, window_size, noverlap, window_size, fs)

end