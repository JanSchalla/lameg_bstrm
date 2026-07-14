function [sFiles, options] = simulate_coherence_MEG2LFP(data_struct, headmodel_fname, sim_params, study_id)
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
%       .sim_loc_extra  – Source Idx of signal A to be coherent
%       .sim_loc_intra  – Source Idx of signal B to be coherent -> sensor idx
%       .woi          – Window of interest [t_start t_end] in seconds.
%                       Either (1×2) for all locations, or (nLocs×2).
%       .fbands          – Frequency of interest [Hz].
%                       Either scalar or (nLocs×1).
%       .DipoleMoment - Dipole strength (nAm) to be simulated. Is converted
%                       to Am. Either a scalar or (nLocs×1).
%       .snr_dB       – Target SNR in dB (applied at sensor level).
%       .targetCoh    - Target Coh (between 0 ant 1).
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
sfreq = 250;
DipoleMoment = 10; 

% Parse inputs
if exist('sim_params', 'var') && ~isempty(sim_params)

    if isfield(sim_params, 'sfreq')
        sfreq = sim_params.sfreq;
    end

    if isfield(sim_params, 'sim_loc_extra')
        sim_loc_extra = sim_params.sim_loc_extra;
    else
        error('No Simulation location on the cortex specified.');
    end

    if isfield(sim_params, 'sim_loc_intra')
        sim_loc_intra = sim_params.sim_loc_intra;
    else
        error('No Simulation location for LFP specified.');
    end

    if isfield(sim_params, 'woi')
        wois = sim_params.woi;
    else
        error('No timewindow to simualte specified.');
    end

    if isfield(sim_params, 'fbands')
        fbands = sim_params.fbands;
    else
        error('No Frequency to simulate specified.');
    end

    if isfield(sim_params, 'DipoleMoment')
        DipoleMoment = sim_params.DipoleMoment;
    end
    
    if isfield(sim_params, 'snr_dB')
        snr_dB = sim_params.snr_dB;
    else
        error('No SNR (dB) specified. Specify via params.snr_dB');
    end

    if isfield(sim_params, 'targetCoh')
        targetCoh = sim_params.targetCoh;
    end
end

% -------------------------------------------------------------------------
%  Initialise outputs
% -------------------------------------------------------------------------
options = struct();
options.sim_params = sim_params;
options.used_headmodel = headmodel_fname;

% -------------------------------------------------------------------------
%  Validate parameter dimensions
% -------------------------------------------------------------------------
n_locs = size(sim_loc_extra, 1);
n_woi = size(wois, 1);
n_fbands = size(fbands, 1);
ns = size(data_struct.Time, 2);

% Convert Dipole moment from nAm to Am
DipoleMoment = DipoleMoment*1e-9;

% Choose target coherence and convert to mixing weight
coh_snr = targetCoh / (1 - targetCoh);  % invert Coh = SNR/(SNR+1)
coh_w = sqrt(coh_snr / (coh_snr + 1));         % weight on shared signal
w_indep = sqrt(1 - coh_w^2);           % weight on independent noise (keeps total power = 1)

if n_locs == 1
    % Single location: both woi and foi must also be scalars / single rows.
    if ~(n_woi == 1 && n_fbands == 1)
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
intra_chans = false(size(grad_chans));
intra_chans(sim_loc_intra) = true; 

head_model = in_bst_headmodel(headmodel_fname);
anatomy = in_tess_bst(head_model.SurfaceFile);
Gain_constrained = bst_gain_orient(head_model.Gain, head_model.GridOrient);

datamat_template = db_template('datamat');
er_template = db_template('datamat');
% -------------------------------------------------------------------------
%  Build zero-padded filename format, e.g. "data_simulation_trial001"
% -------------------------------------------------------------------------
%numDigits = ceil(log10(nTrials + 1));   % +1 avoids log10(0) edge-case
fmt       = sprintf('data_coherence_simulation_trial%%0%dd', 0);

% -------------------------------------------------------------------------
%  Pre-compute per-location signals and store time indices in sim_struct
% -------------------------------------------------------------------------
sim_struct = struct();
dt = 1/sfreq;

% Fixed realizations used both for calibration and final simulation
coh_signal = bandpass(randn(1, ns), fbands, sfreq);
coh_signal = coh_signal / std(coh_signal);

indep_extra = randn(1, ns);   % broadband background, unit variance
indep_lfp   = randn(1, ns);   % broadband background, unit variance

% Solve coh_w such that measured coherence == targetCoh
obj_fun = @(cw) measured_coherence(cw, coh_signal, indep_extra, indep_lfp, fbands, sfreq) - targetCoh;
coh_w = fzero(obj_fun, [1e-4, 0.999]);
w_indep = sqrt(1 - coh_w^2);

% Map woi onto sample indices
[~, sim_struct.t_min_idx] = min(abs(data_struct.Time -min(wois)));
[~, sim_struct.t_max_idx] = min(abs(data_struct.Time -max(wois)));

% ---- Build final signals using the calibrated weight ----
source_ts_unitless = w_indep * indep_extra + coh_w * coh_signal;
LFP_ts = w_indep * indep_lfp + coh_w * coh_signal;

% Assign to sim_struct
sim_struct.signal_extra = zeros(1, ns);
sim_struct.signal_extra(sim_struct.t_min_idx:sim_struct.t_max_idx) = ...
    source_ts_unitless * DipoleMoment;
sim_struct.signal_intra = zeros(1, ns);
sim_struct.signal_intra(sim_struct.t_min_idx:sim_struct.t_max_idx) = ...
    LFP_ts;
sim_struct.foi = fbands;
sim_struct.woi = wois;

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
anatomy.Atlas(user_scout_idx).Scouts(scout_idx).Vertices = sim_loc_extra;
anatomy.Atlas(user_scout_idx).Scouts(scout_idx).Seed = min(sim_loc_extra);
anatomy.Atlas(user_scout_idx).Scouts(scout_idx).Color = [0.2 0.5 0.3];
anatomy.Atlas(user_scout_idx).Scouts(scout_idx).Label = sprintf('Simulation Vertex: %i', sim_loc_extra);

% Save anatomy file
save(file_fullpath(head_model.SurfaceFile), '-struct', 'anatomy');

% -------------------------------------------------------------------------
%  Generate nTrials of simulated sensor data
% -------------------------------------------------------------------------

disp('Starting simulation ...')

% Build coherent signal matrix (n_locs x ns)
coh_signal = vertcat(sim_struct.signal_extra);   % assumes each .signal is 1 x ns

% Compute all active rows at once
%source_ts = w_indep .* sim_struct.noise_extra + coh_w .* coh_signal;  % n_locs x ns

% Single sparse assignment instead of n_locs separate ones
sources = spalloc(size(head_model.GridLoc, 1), ns, n_locs * ns);
sources(sim_loc_extra, :) = sim_struct.signal_extra;

% Forward-project sources to sensor space
sensor = Gain_constrained * sources;

%LFP_ts =  w_indep .* sim_struct.noise_intra + coh_w .* sim_struct.signal_intra;
sensor(sim_loc_intra, :) = sim_struct.signal_intra;

% % -----------------------------------------------------------------
% %  Compute signal power in the woi for SNR-matched noise generation
% %  Use the time window of the first (or only) location as reference.
% % -----------------------------------------------------------------
ref_min = sim_struct(1).t_min_idx;
ref_max = sim_struct(1).t_max_idx;

% -----------------------------------------------------------------
%  Test SPMs way of estimating white noise power to ass
% -----------------------------------------------------------------
% Here i deviate from SPMs apprach and calcualte the rms only over the
% period where a signal is simualted. Otherwise sensor level data is
% shows to be to big by ~1-2 orders of magnitude

std_GRAD = std(sensor(grad_chans, ref_min:ref_max), [], 2); % 1e-11
std_MAG = std(sensor(mag_chans, ref_min:ref_max), [], 2); % 1e-13

rms_GRAD = mean(std_GRAD); % 1e-13
rms_MAG = mean(std_MAG); %1e-14

% Scale noise separately for grad and mag
whitenoise_GRAD = rms_GRAD .* (10^(-snr_dB/20)); %1e-13
whitenoise_MAG = rms_MAG .* (10^(-snr_dB/20)); % 1e-14

noise = zeros(size(sensor));
noise(grad_chans, :) = randn(size(sensor(grad_chans, :))) * whitenoise_GRAD;
noise(mag_chans, :) = randn(size(sensor(mag_chans, :))) * whitenoise_MAG;

% Measure achieved SNR before adding noise
% As SPM calculates noise (on an amplitude basis, this 20*log10())
snr_grad = 20*log10(rms_GRAD/whitenoise_GRAD);
snr_mag = 20*log10(rms_MAG/whitenoise_MAG);

% Add noise; zero out non-MEG channels
sensor = sensor + noise;
sensor(~mag_chans & ~grad_chans & ~intra_chans, :) = 0;

%% Testing Coherence simulation before adding noise on source level
% [cxy_source, F] = mscohere(sim_struct.signal_extra, sim_struct.signal_intra, 256, 128, 256, 250);
% plot(F, cxy_source, 'DisplayName', 'Source Coherence (No Noise)');
% %% Testing Coherence simulation after adding noise on sensor level
% [cxy_sensor, ~] = mscohere(sensor(mag_chans | grad_chans, :)', sensor(intra_chans, :)', hamming(256), 128, 256, 250);
% hold on
% plot(F, mean(cxy_sensor, 2), 'DisplayName', 'Sensor Coherenc (added Noise)');
% legend();
% -----------------------------------------------------------------
%  Populate datamat and save to Brainstorm database
% -----------------------------------------------------------------
foi_str = num2str(sim_struct(1).foi);
woi_str = sprintf('%.3f-%.3f ms', ...
                  min(sim_struct(1).woi) * 1e3, ...
                  max(sim_struct(1).woi) * 1e3);

datamat_template.F            = sensor;
datamat_template.ChannelFlag  = ones(1, size(sensor, 1));
datamat_template.ColormapType = data_struct.ColormapType;
datamat_template.Comment      = sprintf( ...
    'Simulation (%s Hz, Coh: %.2f, SNR: %.2f dB) (#%d)', ...
    foi_str, targetCoh, ...
    mean([snr_grad, snr_mag]), 1);
datamat_template.DataType     = 'recordings';
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


fname = sprintf(fmt, 1);
fname_full = fullfile(cond_path, fname);
sFiles = file_short(fname_full);
    
save(fname_full, '-struct', "datamat_template");
db_add_data(study_id, sFiles, datamat_template);

%% save empty room
er_template.F = noise;
er_template.ChannelFlag  = ones(1, size(sensor, 1));
er_template.ColormapType = data_struct.ColormapType;
er_template.Comment      = 'Simulated Empty Room';
er_template.DataType     = 'recordings';
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
er_template.Time         = data_struct.Time;

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

%% helper functions
function c = measured_coherence(coh_w, coh_signal, indep_a, indep_b, fband, sfreq)
    w_indep = sqrt(1 - coh_w^2);
    x = w_indep*indep_a + coh_w*coh_signal;
    y = w_indep*indep_b + coh_w*coh_signal;
    [Cxy, F] = mscohere(x, y, hamming(256), 128, 256, sfreq);
    band_mask = F >= fband(1) & F <= fband(2);
    c = mean(Cxy(band_mask));
end
end