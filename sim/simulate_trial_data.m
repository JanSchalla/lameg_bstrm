function [sFiles, options] = simulate_trial_data(data_struct, headmodel_fname, sim_params, study_id)
% SIMULATE_TRIAL_DATA  Simulate MEG trial data and save to a Brainstorm study.
%
%   [sFiles, options] = simulate_trial_data(data_struct, headmodel_fname,
%                                           sim_params, study_id)
%
%   Generates nTrials of synthetic MEG data by projecting sinusoidal source
%   signals through a forward model and adding channel-type-specific
%   Gaussian white noise at a requested SNR.  Results are written to the
%   Brainstorm database and the corresponding file list is returned.
%
% INPUTS
%   data_struct       – Brainstorm datamat used as a template (fields: Time,
%                       ColormapType, DataType, Device, DisplayUnits, Leff,
%                       nAvg, Std).
%   headmodel_fname   – Short (relative) Brainstorm path to the head-model
%                       file that provides Gain and GridOrient / GridLoc.
%   sim_params        – Struct with simulation settings:
%       .nTrials      – Number of trials to generate.
%       .sfreq        – Sampling frequency [Hz].
%       .sim_loc      – (nLocs × nVerts) matrix of source vertex indices.
%       .woi          – Window of interest [t_start t_end] in seconds.
%                       Either (1×2) for all locations, or (nLocs×2).
%       .foi          – Frequency of interest [Hz].
%                       Either scalar or (nLocs×1).
%       .DipoleMoment - Dipole moment (nAm) to be simulated. Is converted
%                       to Am. Either a scalar or (nLocs×1).
%       .snr_dB       – Target SNR in dB (applied at sensor level).
%   study_id          – Brainstorm study index (integer).
%
% OUTPUTS
%   sFiles            – (1 × nTrials) cell array of short Brainstorm file
%                       paths for the saved trials.
%   options           – Struct echoing the simulation configuration plus
%                       per-trial SNR diagnostics:
%       .sim_params       – Copy of input sim_params.
%       .used_headmodel   – headmodel_fname.
%       .snr_grad         – (1 × nTrials) achieved gradiometer SNR [dB].
%       .snr_mag          – (1 × nTrials) achieved magnetometer SNR [dB].
%
% NOTES
%   Source amplitudes are randomised trial-by-trial with a Gaussian offset
%   (mean ≈ 1.87 nAm, std ≈ 1.16 nAm).  Only MEG MAG and MEG GRAD channels
%   carry signal; all other channels are set to zero.
%
% EXAMPLE
%   sim_params.nTrials = 50;
%   sim_params.sfreq   = 1000;
%   sim_params.sim_loc = [42 43 44];
%   sim_params.woi     = [0.1 0.4];
%   sim_params.foi     = 20;
%   sim_params.snr_dB  = 0;
%   [sFiles, opts] = simulate_trial_data(data_struct, hm_fname, sim_params, 1);

% -------------------------------------------------------------------------
%  Defaults
% -------------------------------------------------------------------------

nTrials = 540;
sfreq = 250;
DipoleMoment = 10; 

% Parse inputs
if exist('sim_params', 'var') && ~isempty(sim_params)
    if isfield(sim_params, 'nTrials')
        nTrials = sim_params.nTrials;
    end

    if isfield(sim_params, 'sfreq')
        sfreq = sim_params.sfreq;
    end

    if isfield(sim_params, 'sim_loc')
        sim_loc = sim_params.sim_loc;
    else
        error('No Simulation location specified.');
    end

    if isfield(sim_params, 'woi')
        wois = sim_params.woi;
    else
        error('No timewindow to simualte specified.');
    end

    if isfield(sim_params, 'foi')
        fois = sim_params.foi;
    else
        error('No Frequency to simulate specified.');
    end

    if isfield(sim_params, 'DipoleMoment')
        DipoleMoment = sim_params.DipoleMoment;
    end
    
    if isfield(sim_params, 'snr_dB')
        snr_dB = sim_params.snr_dB;
    else
        error('No SNR (dB) specified.');
    end
end

% -------------------------------------------------------------------------
%  Initialise outputs
% -------------------------------------------------------------------------
sFiles = cell(1, nTrials);
options = struct();
options.sim_params = sim_params;
options.used_headmodel = headmodel_fname;

% -------------------------------------------------------------------------
%  Validate parameter dimensions
% -------------------------------------------------------------------------
n_locs = size(sim_loc, 1);
n_woi = size(wois, 1);
n_foi = size(fois, 1);
ns = size(data_struct.Time, 2);

% Convert Dipole moment from nAm to Am
DipoleMoment = DipoleMoment*1e-9;

empty_room = zeros(size(data_struct.F, 1), ns*floor(nTrials/10));

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

% -------------------------------------------------------------------------
%  Load Brainstorm study / channel / head-model data
% -------------------------------------------------------------------------
sStudy = bst_get('Study', study_id);
cond_path = fileparts(file_fullpath(sStudy.FileName));
channels = in_bst_channel(sStudy.Channel.FileName);

mag_chans = strcmp({channels.Channel.Type}, 'MEG MAG')';
grad_chans = strcmp({channels.Channel.Type}, 'MEG GRAD')';

head_model = in_bst_headmodel(headmodel_fname);
anatomy = in_tess_bst(head_model.SurfaceFile);
Gain_constrained = bst_gain_orient(head_model.Gain, head_model.GridOrient);

datamat_template = db_template('datamat');
er_template = db_template('datamat');
% -------------------------------------------------------------------------
%  Build zero-padded filename format, e.g. "data_simulation_trial001"
% -------------------------------------------------------------------------
numDigits = ceil(log10(nTrials + 1));   % +1 avoids log10(0) edge-case
fmt       = sprintf('data_simulation_trial%%0%dd', numDigits);

% -------------------------------------------------------------------------
%  Pre-compute per-location signals and store time indices in sim_struct
% -------------------------------------------------------------------------
sim_struct = struct();
dt = 1/sfreq;

for i = 1:n_locs
    
    if n_woi ~= 1
        woi = wois(i, :);
    else
        woi = wois;
    end

    % Map woi onto sample indices
    [~, sim_struct(i).t_min_idx] = min(abs(data_struct.Time -min(woi)));
    [~, sim_struct(i).t_max_idx] = min(abs(data_struct.Time -max(woi)));

    % Build sinusoidal signal for the woi
    time_diff = diff(woi);
    t = (0:dt:time_diff);
    
    if numel(fois) ~= 1
        foi = fois(i);
    else
        foi = fois;
    end

    % Reconcile rounding differences between time vector and sample window
    woi_samples = sim_struct(i).t_min_idx : sim_struct(i).t_max_idx;
    if length(woi_samples) < length(t)
        warning(['simulate_trial_data: simulated signal (loc %d) is longer ' ...
                 'than the sample window — signal will be cropped.'], i);
        t = t(1 : length(woi_samples));
    elseif length(woi_samples) > length(t)
        warning(['simulate_trial_data: sample window (loc %d) is longer '  ...
                 'than the simulated signal — window end will be cropped.'], i);
        sim_struct(i).t_max_idx = t_min_idx + length(t) - 1;
    end

    sim_struct(i).signal = zeros(1, ns);
    sim_struct(i).signal(sim_struct(i).t_min_idx : sim_struct(i).t_max_idx) = ...
        sin(2 * pi * foi * t) * DipoleMoment;
    sim_struct(i).foi = foi;
    sim_struct(i).woi = woi;

    % Insert sim locations as scouts to anatomy
    user_scout_idx = find(ismember({anatomy.Atlas.Name}, 'User scouts'));
    if isempty(anatomy.Atlas(user_scout_idx).Scouts)
        scout_idx = 1;
    else
        if isempty(anatomy.Atlas(user_scout_idx).Scouts(1).Vertices)
            scout_idx = 1;
        else
            scout_idx = size(anatomy.Atlas(user_scout_idx).Scouts, 2) + 1; 
        end
    end
    anatomy.Atlas(user_scout_idx).Scouts(scout_idx) = db_template('Scout');
    anatomy.Atlas(user_scout_idx).Scouts(scout_idx).Vertices = sim_loc(i);
    anatomy.Atlas(user_scout_idx).Scouts(scout_idx).Seed = min(sim_loc(i));
    anatomy.Atlas(user_scout_idx).Scouts(scout_idx).Color = [0.2 0.5 0.3];
    anatomy.Atlas(user_scout_idx).Scouts(scout_idx).Label = sprintf('Simulation Vertex: %i', sim_loc(i));
end

% Save anatomy file
save(file_fullpath(head_model.SurfaceFile), '-struct', 'anatomy');

% -------------------------------------------------------------------------
%  Generate nTrials of simulated sensor data
% -------------------------------------------------------------------------
snr_grad = zeros(1, nTrials);
snr_mag = zeros(1, nTrials);

disp('Starting simulation ...')
fprintf('Progress: %3d%%\n', 0);
start_idx = 1;

for iTrial = 1:nTrials

    fprintf(1, '\b\b\b\b%3.0f%%', 100*iTrial/nTrials);
    
    % Initialise source matrix (nSources × nSamples)
    sources = sparse(size(head_model.GridLoc, 1), length(data_struct.Time));
    
     % Assign (randomly scaled) signal to each active source
    for iLoc = 1 : n_locs
        %src_scaling = 1.1616e-9 + 1.8735e-9; 
        %src_scaling = randn(length(sim_params.sim_loc(iLoc, :)), 1) ...
        %              * 1.1616e-9 + 1.8735e-9; % Parameter eyeballed to achieve realisitic sensor level scaling after projecting
        sources(sim_loc(iLoc, :), :) = sim_struct(iLoc).signal; 
    end
        
    % Forward-project sources to sensor space
    sensor = sparse(Gain_constrained * sources);
    
    % % -----------------------------------------------------------------
    % %  Compute signal power in the woi for SNR-matched noise generation
    % %  Use the time window of the first (or only) location as reference.
    % % -----------------------------------------------------------------
    % ref_min = sim_struct(1).t_min_idx;
    % ref_max = sim_struct(1).t_max_idx;
    % signal_power_grad = mean(mean(abs(sensor(grad_chans, ref_min:ref_max)), 'omitmissing'));
    % signal_power_mag = mean(mean(abs(sensor(mag_chans, ref_min:ref_max)), 'omitmissing'));
    % 
    % % Convert target SNR from dB to linear scale and derive noise std
    % linear_snr = 10^(snr_dB/10);
    % noise_std_grad = signal_power_grad/linear_snr;
    % noise_std_mag = signal_power_mag/linear_snr;
    % 
    % % Draw channel-type-specific Gaussian white noise
    % noise = randn(size(sensor, 1), size(sensor, 2));
    % noise_power = mean(mean(abs(noise)));
    % noise(grad_chans, :) = noise_std_grad * noise(grad_chans, :) / noise_power;
    % noise(mag_chans, :) = noise_std_mag * noise(mag_chans, :) / noise_power;

    % -----------------------------------------------------------------
    %  Test SPMs way of estimating white noise power to ass
    % -----------------------------------------------------------------
    % Here i deviate from SPMs apprach and calcualte the rms only over the
    % period where a signal is simualted. Otherwise sensor level data is
    % shows to be to big by ~1-2 orders of magnitude
    std_GRAD = std(sensor(grad_chans, ref_min:ref_max), [], 2);
    std_MAG = std(sensor(mag_chans, ref_min:ref_max), [], 2);

    rms_GRAD = mean(std_GRAD);
    rms_MAG = mean(std_MAG);

    % Scale noise separately for grad and mag
    whitenoise_GRAD = rms_GRAD .* (10^(-snr_dB/20));
    whitenoise_MAG = rms_MAG .* (10^(-snr_dB/20));

    noise = zeros(size(sensor));
    noise(grad_chans, :) = randn(size(sensor(grad_chans, :))) * whitenoise_GRAD;
    noise(mag_chans, :) = randn(size(sensor(mag_chans, :))) * whitenoise_MAG;

    if mod(iTrial, 10) == 0
        empty_room(:, start_idx:start_idx + ns -1) = noise;
        start_idx = start_idx + ns;
    end
    
    % Measure achieved SNR before adding noise
    % snr_grad(iTrial) = 10*log10(signal_power_grad/(mean(mean(abs(noise(grad_chans, ref_min:ref_max))))));
    % snr_mag(iTrial) = 10*log10(signal_power_mag/(mean(mean(abs(noise(mag_chans, ref_min:ref_max))))));
    % As SPM calculates noise (on an amplitude basis, this 20*log10())
    snr_grad(iTrial) = 20*log10(rms_GRAD/whitenoise_GRAD);
    snr_mag(iTrial) = 20*log10(rms_MAG/whitenoise_MAG);

    % Add noise; zero out non-MEG channels
    sensor = full(sensor) + noise;
    sensor(~mag_chans & ~grad_chans, :) = 0;
    
    % -----------------------------------------------------------------
    %  Populate datamat and save to Brainstorm database
    % -----------------------------------------------------------------
    % Build a human-readable comment (safe for both scalar and vector params)
    foi_str = num2str(sim_struct(1).foi);
    woi_str = sprintf('%.3f-%.3f ms', ...
                      min(sim_struct(1).woi) * 1e3, ...
                      max(sim_struct(1).woi) * 1e3);

    datamat_template.F            = sensor;
    datamat_template.ChannelFlag  = ones(1, size(sensor, 1));
    datamat_template.ColormapType = data_struct.ColormapType;
    datamat_template.Comment      = sprintf( ...
        'Simulation (%s Hz, WOI: %s, SNR: %.2f dB) (#%d)', ...
        foi_str, woi_str, ...
        mean([snr_grad(iTrial), snr_mag(iTrial)]), iTrial);
    datamat_template.DataType     = data_struct.DataType;
    datamat_template.Device       = data_struct.Device;
    datamat_template.DisplayUnits = data_struct.DisplayUnits;
    datamat_template.Events       = db_template('event');
    datamat_template.History      = {'simulate', ...
                                     char(datetime('now', 'Format', ...
                                         'yyyy-MM-dd''T''HH:mm:ss')), ...
                                     datamat_template.Comment};
    datamat_template.Leff         = data_struct.Leff;
    datamat_template.nAvg         = data_struct.nAvg;
    datamat_template.Std          = data_struct.Std;
    datamat_template.Time         = data_struct.Time;


    fname = sprintf(fmt, iTrial);
    fname_full = fullfile(cond_path, fname);
    sFiles{iTrial} = file_short(fname_full);
        
    save(fname_full, '-struct', "datamat_template");
    db_add_data(study_id, sFiles{iTrial}, datamat_template);
end % iTrial

% save empty room
er_template.F = empty_room;
er_template.ChannelFlag  = ones(1, size(sensor, 1));
er_template.ColormapType = data_struct.ColormapType;
er_template.Comment      = 'Simulated Empty Room';
er_template.DataType     = data_struct.DataType;
er_template.Device       = data_struct.Device;
er_template.DisplayUnits = data_struct.DisplayUnits;
er_template.Events       = db_template('event');
er_template.History      = {'simulate', ...
                                 char(datetime('now', 'Format', ...
                                     'yyyy-MM-dd''T''HH:mm:ss')), ...
                                 datamat_template.Comment};
er_template.Leff         = data_struct.Leff;
er_template.nAvg         = data_struct.nAvg;
er_template.Std          = data_struct.Std;
er_template.Time         = (0:size(empty_room, 2)-1) / sfreq;

fname_full = fullfile(cond_path, 'data_empty_room.mat');
save(fname_full, '-struct', "er_template");
%db_add_data(study_id, fname_full);

% Reload the study so the new trials appear in the Brainstorm GUI
db_reload_studies(study_id, 1);

% -------------------------------------------------------------------------
%  Store diagnostics and print summary
% -------------------------------------------------------------------------
options.snr_grad = snr_grad;
options.snr_mag = snr_mag;

fprintf('\nSimulation complete!\n');
fprintf('Signal simulated with a mean SNR of %.2f dB (target: %.2f dB)\n', mean([mean(snr_grad), mean(snr_mag)]), snr_dB);
fprintf('Mean Gradiometer SNR: %.2f dB\n', mean(snr_grad));
fprintf('Mean Magnetometer SNR: %.2f dB\n', mean(snr_mag));

