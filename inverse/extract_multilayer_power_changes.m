function extract_multilayer_power_changes(sFiles, protocol_name, params)

if ~brainstorm('status')
    brainstorm nogui
end

% Parse inputs
if exist('params', 'var') && ~isempty(params)
    if isfield(params, 'prompt')
        answer = params.prompt;
    end
end

 % Open protocol
protocol = bst_get('Protocol', protocol_name);
gui_brainstorm('SetCurrentProtocol', protocol);

tf_template = db_template('timefreq');

% sFiles cell of relaitves pathes to the trial data from the brainstorm db
sData = in_bst_data(sFiles{1});
time = sData.Time;
clear sData
%% 1st: Let user define baseline, window of interest (woi) and frequency of interest (foi).
if ~exist("answer", "var")
    dlgTitle = 'Specify multilayer extraction parameter';
    fieldsize = [10 50];
    prompt = 'freq id/ freq range / BL win / WOI';
    example = {'theta / 4, 8 / -0.5, -0.1 / 0, 1'};
    answer = inputdlg(prompt, dlgTitle, fieldsize, example);

    jobs = cellstr(answer{:});
else
    jobs = cellstr(answer);
end


sfreq = round(1/(time(2) - time(1)));
%% Check if all files are of the same run
nTrialsTotal = length(sFiles);
ses_ids = unique(cellfun(@fileparts, sFiles, 'UniformOutput', false));
% group files in sessions
session_files = cell(1, length(ses_ids));
for ses = 1:length(ses_ids)
    file_mask = cellfun(@(x) contains(x, ses_ids{ses}), sFiles);
    session_files{ses} = {sFiles{file_mask}};
end

fprintf('%i Sessions identfied.\n', length(ses_ids));

for ses=1:length(ses_ids)
    fprintf('Working on Session: %i\n', ses);

    % init brainstorm
    bstrm_path = session_files{ses}{1};
    if length(split(bstrm_path, "/")) > length(split(bstrm_path, "\"))
        separator = "/";
    else
        separator = "\";
    end
    tokens = split(bstrm_path, separator);
    protocol_id = find(ismember(tokens, "brainstorm_db"), 1) + 1;
    % check if previous path was relative
    if isempty(protocol_id)
        bstrm_path = file_fullpath(bstrm_path);
        % Redo if previous path was relative
        tokens = split(bstrm_path, separator);
    end

    subject_id = find(ismember(tokens, "data"), 1) + 1;
    subject_name = convertStringsToChars(tokens{subject_id});
    
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
    
    % Populate output structure
    ses_template = tf_template;
    ses_template.HeadModelFile = sKernel.HeadModelFile;
    ses_template.SurfaceFile = sKernel.SurfaceFile;
    ses_template.RowNames = [1:nSources];
    ses_template.Comment = ses_ids{ses};
    ses_template.TF = zeros(nSources, 1, length(jobs));
    ses_template.Time = [min(time), max(time)];
    ses_template.Measure = 'Power';
    ses_template.Method = 'hilbert';
    ses_template.nAvg = 1;
    ses_template.Leff = 1;

    disp('Starting extraction ...')
    fprintf('Progress: %3d%%\n', 0);
    for iTrial = 1:sesTrials
        fprintf(1, '\b\b\b\b%3.0f%%', 100*(iTrial/sesTrials));

        trial_template = ses_template;

        % Load data for each trial
        [~, trial_id] = fileparts(session_files{ses}{iTrial});
        sRawData = in_bst_data(session_files{ses}{iTrial});
    
        % Reduce sensor data to only used channels in source reconstruction
        sensorData = sRawData.F(sKernel.GoodChannel, :);

        for job=1:length(jobs)
        
            %% Initialize current extraction params
            job_specs = split(jobs{job}, '/');
            job_specs = cellfun(@(x) x(find(~isspace(x))), job_specs, 'UniformOutput', false);
        
            job_struct = struct( ...
                'id', job_specs{1}, ...
                'freq_range', str2double(split(job_specs{2}, ',')), ...
                'base_win', str2double(split(job_specs{3}, ',')), ...
                'woi', str2double(split(job_specs{4}, ',')));
        
            if min(job_struct.base_win) < min(time)
                error('Baseline starts outside of trial. Baseline cannot be defined smaller then %.3f seconds (Job ID: %s)', min(time), job_struct.id);
            end
        
            if max(job_struct.base_win) > min(job_struct.woi)
                warning('Baseline and Window if Interest overlap. This is not recommended (Job ID: %s).', job_struct.id);
                dlgTitle    = 'Overlap of Baseline and WOI';
                dlgQuestion = 'Do you wish to continue?';
                choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'No');
            
                if strcmp(choice, 'No')
                    error('Baseline and WOI overlap.');
                end
            end
        
            if max(job_struct.woi) > max(time)
                error('Window of Interest ends outside of trial bound. WOI cannot be defined bigger then %.3f seconds (Job ID: %s)', max(time), job_struct.id);
            end
            
            if max(job_struct.freq_range) > sfreq/2
                error('Selected Frequency of Interest is bigger then Nyquist Frequency (%i Hz).', round(sfreq/2));
            end
                        
            if ~isfield(trial_template, 'Freqs')
                % Update Frequency Bins
                trial_template.Freqs{1, 1} = char(job_struct.id);
                trial_template.Freqs{1, 2} = char(join(string(job_struct.freq_range), ','));
                trial_template.Freqs{1, 3} = char(sprintf('BL: %.2f-%.2f; WOI: %.2f-%.2f', job_struct.base_win, job_struct.woi));
                
                % Update history
                trial_template.History{1, 1} = char(datestr(now, 'dd/mm/yy-HH:MM'));
                trial_template.History{1, 2} = char('compute'); 
                trial_template.History{1, 3} = char(sprintf('extract_multilayer_power_changes | %s; BL: %f %f; WOI: %f %f; Freq: %i %i', job_struct.id, job_struct.base_win, job_struct.woi, job_struct.freq_range));
            else 
                freq_id = size(trial_template.Freqs, 1);
                % Update Frequency bins
                trial_template.Freqs{freq_id+1, 1} = char(job_struct.id);
                trial_template.Freqs{freq_id+1, 2} = char(join(string(job_struct.freq_range), ','));
                trial_template.Freqs{freq_id+1, 3} = char(sprintf('BL: %.2f-%.2f; WOI: %.2f-%.2f', job_struct.base_win, job_struct.woi));
                
                % Update history
                trial_template.History{freq_id+1, 1} = char(datestr(now, 'dd/mm/yy-HH:MM'));
                trial_template.History{freq_id+1, 2} = char('compute'); 
                trial_template.History{freq_id+1, 3} = char(sprintf('extract_multilayer_power_changes | %s; BL: %f %f; WOI: %f %f; Freq: %i %i', job_struct.id, job_struct.base_win, job_struct.woi, job_struct.freq_range));
            end
            
            % Bandpasss filter sensor data
            sensorData = bandpass(sensorData', job_struct.freq_range, sfreq)';
            
            % bring time information in sample space
            [~, base_min_idx] = min(abs(time - job_struct.base_win(1)));
            [~, base_max_idx] = min(abs(time - job_struct.base_win(2)));
            base_win_samples = [base_min_idx base_max_idx];
            clear('base_min_idx', "base_max_idx");
            
            [~, woi_min_idx] = min(abs(time - job_struct.woi(1)));
            [~, woi_max_idx] = min(abs(time - job_struct.woi(2)));
            woi_samples = [woi_min_idx woi_max_idx];
            clear('woi_min_idx', "woi_max_idx");

            % extract job specific data
            source_data_base = sKernel.ImagingKernel * sensorData(:, base_win_samples(1):base_win_samples(2));
            source_data_woi = sKernel.ImagingKernel * sensorData(:, woi_samples(1):woi_samples(2));
           
            % In the original lameg code (from Jimmy) the absolute hilbert
            % is not squared. Here we square to have an assessment of
            % instantaneous power.
            source_power_base = mean(abs(hilbert(source_data_base)).^2, 2);
            source_power_woi = mean(abs(hilbert(source_data_woi)).^2, 2);
            
            power_change = source_power_woi - source_power_base;
            
            trial_template.TF(:, 1, job) = power_change;

            trial_template.Comment = char(join([trial_template.Comment trial_id "multilayer power change"], " | "));
            trial_template.History{end+1, 1} = char(datestr(now, 'dd/mm/yy-HH:MM'));
            trial_template.History{end, 2} = char('compute'); 
            trial_template.History{end, 3} = char(sprintf('extract_multilayer_power_changes | %s; BL: %f %f; WOI: %f %f; Freq: %i %i', job_struct.id, job_struct.base_win, job_struct.woi, job_struct.freq_range));
        end 
        
        % Link it to the trial sensor data and to the multilayer kernel!
        trial_template.DataFile = char(sprintf('link|%s|%s', file_short(char(kernel_file)), session_files{ses}{iTrial}));
        
        out_fname = sprintf('timefreq_%s_psd_multilayer_power_change.mat', trial_id(end-7:end));
        
        save(fullfile(bstrm_out_path, out_fname), '-struct', 'trial_template');
    end
end

fpintf('\nExtraction finished!\n');