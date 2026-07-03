function freq_csd = bst_fT_calculateCSDMatrix(ftData, params)

%% ------------------------------------------------------------------
%% Defaults
%% ------------------------------------------------------------------
freq_id = 'None';
freq_range = [];
data_cov = false;
noise_cov = false;

if exist('params', 'var') && ~isempty(params)
    if isfield(params, 'freq_id')
        freq_id = params.freq_id;
    end

    if isfield(params, 'freq_range')
        freq_range = params.freq_range;
    end

    if isfield(params, 'data_cov')
        data_cov = params.data_cov;
    end
    
    if isfield(params, 'noise_cov')
        noise_cov = params.noise_cov;
    end
end

% Verify input
if isempty(freq_range)
    error('No Frequency range to calculate cross spectral density supplied (freq_range).');
end

if strcmp(freq_id, 'None')
    warning('No Frequency Identifier (freq_id) specified. Data will be non-identifiable.')
end

if data_cov
    warning('Fieldtrips CFG will be saved as Data Covariance');
elseif noise_cov
    warning('Fieldtrips CFG will be saved as Noise Covariance');
else
    error('Not specified, if data covariance (data_cov) or noise covariance (noise_cov) is calculated.');
end

%% ------------------------------------------------------------------
%% Step 3: Compute cross-spectral density matrix per frequency band
%% ------------------------------------------------------------------

centerFreq = mean(freq_range);
halfBW = diff(freq_range)/2;

cfg = [];
cfg.output    = 'powandcsd';
cfg.method    = 'mtmfft';
cfg.taper     = 'hanning';   % single taper - segments already provide averaging (see note above)
cfg.foi       = centerFreq;
cfg.tapsmofrq = halfBW;      % spectral smoothing matched to band half-width
cfg.pad         = 'nextpow2'; % From J. Hirschmanns script
cfg.keeptrials = 'no';       % average CSD across all trials directly
cfg.channel    = ftData.label;          % all MEG + LFP channels
cfg.channelcmb = {ftData.label, ftData.label};  % all pairwise combinations

freq_csd = ft_freqanalysis(cfg, ftData);

fprintf('Computed CSD for band %s (center %.1f Hz, smoothing +/-%.1f Hz)\n', ...
    freq_id, centerFreq, halfBW);

% % Save multilayer results
% ResultsMat = db_template('noisecovmat');
% ResultsMat.NoiseCov = freq_csd;
% 
% if data_cov
%     OutputFile = fullfile(fileparts(DataFile), ...
%         sprintf('ndatacov_CSD_fieldtrip_%s.mat', freq_id));
%     covariance_comment = sprintf('Data Cross Spectral Density (%s, Fieldtrips cfg)', freq_id);
% elseif noise_cov
%     OutputFile = fullfile(fileparts(DataFile), ...
%         sprintf('noisecov_CSD_fieldtrip_%s.mat', freq_id));
%     covariance_comment = sprintf('Noise Cross Spectral Density (%s, Fieldtrips cfg)', freq_id);
% end
% 
% ResultsMat.Comment = covariance_comment;
% 
% bst_save(OutputFile, ResultsMat, 'v6');
% 
% fprintf('Saved %s FieldTrip CSD to: %s\n', freq_id, OutputFile);
