function [sFiles, options] = simulate_trial_data(data_struct, headmodel_fname, sim_params, study_id)

sFiles = cell(1, sim_params.nTrials);
options = struct();
options.sim_params = sim_params;
options.used_headmodel = headmodel_fname;

% Check sim_params to be congruent
n_locs = size(sim_params.sim_loc, 1);
n_woi = size(sim_params.woi, 1);
n_foi = size(sim_params.foi, 1);
ns = size(data_struct.Time, 2);

if n_locs == 1
    % Single location: both woi and foi must also be scalars / single rows.
    if ~(n_woi == 1 && n_foi == 1)
        error(['simulate_trial_data: only one location specified, ' ...
               'but multiple woi/foi entries were found.']);
    end
else
    % Multiple locations: either one shared param set, or one per location.
    if ~((n_woi == 1 && n_foi == 1) || (n_woi == n_locs && n_foi == n_locs))
        error(['simulate_trial_data: for %d locations, woi and foi must ' ...
               'each be either (1×…) or (%d×…).'], n_locs, n_locs);
    end
end

% Load 
sStudy = bst_get('Study', study_id);
cond_path = fileparts(file_fullpath(sStudy.FileName));
channels = in_bst_channel(sStudy.Channel.FileName);
mag_chans = strcmp({channels.Channel.Type}, 'MEG MAG')';
grad_chans = strcmp({channels.Channel.Type}, 'MEG GRAD')';

head_model = in_bst_headmodel(headmodel_fname);
Gain_constrained = bst_gain_orient(head_model.Gain, head_model.GridOrient);

%anatomy = load(file_fullpath(head_model.SurfaceFile));
datamat_template = db_template('datamat');

% Create filename with leading zeros, depending on how much trials are
% simulated
numDigits = ceil(log10(sim_params.nTrials));
fmt = sprintf('data_simulation_trial%%0%dd', numDigits);

% Get WOI onset and offset in trial samples
sim_struct = struct();
for i = 1:n_locs
    
    if n_woi ~= 1
        woi = sim_params.woi(i, :);
    else
        woi = sim_params.woi;
    end

    [~, t_min_idx] = min(abs(data_struct.Time -min(woi)));
    [~, t_max_idx] = min(abs(data_struct.Time -max(woi)));

    sim_struct(i).t_min_idx = t_min_idx;
    sim_struct(i).t_max_idx = t_max_idx;

    % Simulate the signal
    dt = 1/sim_params.sfreq;
    time_diff = diff(woi);
    t = (0:dt:time_diff);
    
    if numel(sim_params.foi) ~= 1
        foi = sim_params.foi(i);
    else
        foi = sim_params.foi;
    end

    % Check for mismatch in length of signal
    if length(t_min_idx:t_max_idx) < length(t)
        warning('Length of simulated signal is longer then the time specification. Signal will be cropped ...');
        t = t(1:length(t_min_idx:t_max_idx));
    elseif length(t_min_idx:t_max_idx) > length(t)
        warning('Length of simulated signal is shorter then the time specification. Time specification will be cropped ...');
        t_max_idx = t_max_idx - (length(t_min_idx:t_max_idx) - length(t));
    end

    sim_struct(i).signal = zeros(1, ns);
    sim_struct(i).signal(t_min_idx:t_max_idx) = sin(2*pi*foi*t);
end

% Create nTrials of simulated data
snr_grad = zeros(1, sim_params.nTrials);
snr_mag = zeros(1, sim_params.nTrials);

disp('Starting simulation ...')
fprintf('Progress: %3d%%\n', 0);
for iTrial = 1:sim_params.nTrials
    % Create and progress update
    fprintf(1, '\b\b\b\b%3.0f%%', 100*iTrial/sim_params.nTrials);
    
    % Preproc
    sources = zeros(size(head_model.GridLoc, 1), length(data_struct.Time));
    
    % Assign signal to the sources
    for iLoc = 1:n_locs
        src_scaling = randn(length(sim_params.sim_loc(iLoc, :)), 1) * 1.1616e-9 + 1.8735e-9;
        sources(sim_params.sim_loc(iLoc, :), :)= sim_struct(iLoc).signal .* src_scaling;
    end
        
    % Project sources to sensors
    sensor = Gain_constrained * sources;
    
    % In Bonaiutos Simulation this is only done at the sensor level
    % Create Noise
    % Separate for each channeltype
    ref_min = sim_struct(1).t_min_idx;
    ref_max = sim_struct(1).t_max_idx;
    signal_power_grad = mean(mean(sensor(grad_chans, ref_min:ref_max).^2, 'omitmissing'));
    signal_power_mag = mean(mean(sensor(mag_chans, ref_min:ref_max).^2, 'omitmissing'));
    
    % Transfrom SNR in dB to a linear scale
    linear_snr = 10^(sim_params.snr_dB/10);
    noise_std_grad = sqrt(signal_power_grad/linear_snr);
    noise_std_mag = sqrt(signal_power_mag/linear_snr);

    % Gaussian white noise
    noise = randn(size(sensor, 1), size(sensor, 2));
    noise(grad_chans, :) = noise_std_grad * noise(grad_chans, :);
    noise(mag_chans, :) = noise_std_mag * noise(mag_chans, :);
    
    % Calculate SNR on the vertex of interest
    snr_grad(iTrial) = 10*log10(signal_power_grad/(mean(mean(noise(grad_chans, min_time_idx:max_time_idx).^2))));
    snr_mag(iTrial) = 10*log10(signal_power_mag/(mean(mean(noise(mag_chans, min_time_idx:max_time_idx).^2))));

    % Add noise to the signal
    sensor = sensor + noise;

    % Set everything which is not MEG to 0.
    sensor(~mag_chans & ~grad_chans, :) = 0;
    
    %plot(data_struct.Time, sensor)
    
    datamat_template.F = sensor;
    datamat_template.ChannelFlag = ones(1, size(sensor, 1));
    datamat_template.ColormapType = data_struct.ColormapType;
    datamat_template.Comment = sprintf('Simulation (%iHz, WOI: %.3f-%.3fms, SNR: %.2f dB) (#%i)', sim_params.foi, min(sim_params.woi), max(sim_params.woi), mean([snr_grad(iTrial), snr_mag(iTrial)]), iTrial);
    datamat_template.DataType = data_struct.DataType;
    datamat_template.Device = data_struct.Device;
    datamat_template.DisplayUnits = data_struct.DisplayUnits;
    datamat_template.Events = db_template('event');
    datamat_template.History = {convertStringsToChars(string(datetime('now'))), 'simulate'};
    datamat_template.Leff = data_struct.Leff;
    datamat_template.nAvg = data_struct.nAvg;
    datamat_template.Std = data_struct.Std;
    datamat_template.Time = data_struct.Time;

    fname = sprintf(fmt, iTrial);
    fname_full = fullfile(cond_path, fname);
    sFiles{iTrial} = file_short(fname_full);
        
    save(fname_full, '-struct', "datamat_template");
    db_add_data(study_id, sFiles{iTrial}, datamat_template);
end
% Reload study to make trials appear
db_reload_studies(study_id, 1);

options.snr_grad = snr_grad;
options.snr_mag = snr_mag;

fprintf('\nSimulation complete!\n');
fprintf('Signal simulated with a mean SNR of %.2f dB (target: %.2f dB)\n', mean([mean(snr_grad), mean(snr_mag)]), sim_params.snr_dB);
fprintf('Mean Gradiometer SNR: %.2f dB\n', mean(snr_grad));
fprintf('Mean Magnetometer SNR: %.2f dB\n', mean(snr_mag));

