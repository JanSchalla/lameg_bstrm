function freq_csd = bst_fT_calculateCSDMatrix(DataFile, extracranialIdx, intracranialIdx, params)

%% ------------------------------------------------------------------
%% Defaults
%% ------------------------------------------------------------------

seg_length_sec = 2;
seg_overlap = 0.5;
freq_id = 'None';
freq_range = [];
data_cov = false;
noise_cov = false;


if exist('params', 'var') && ~isempty(params)
    if isfield(params, 'seg_length_sec')
        seg_length_sec = params.seg_length_sec;
    end

    if isfield(params, 'seg_overlap')
        seg_overlap = params.seg_overlap;
    end

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

if seg_overlap < 0 || seg_overlap > 0.99
    error('Window overlap is not possible. Specify between 0 - 0.99 (seg_overlap).');
end

if data_cov
    warning('Fieldtrips CFG will be saved as Data Covariance');
elseif noise_cov
    warning('Fieldtrips CFG will be saved as Noise Covariance');
else
    error('Not specified, if data covariance (data_cov) or noise covariance (noise_cov) is calculated.');
end

%% ------------------------------------------------------------------
%% Step 0: Initiate toolboxes
%% ------------------------------------------------------------------
%% Load brainstorm project and file

% Parse path
if length(split(DataFile, "/")) > length(split(DataFile, "\"))
    separator = "/";
else
    separator = "\";
end

tokens = split(DataFile, separator);

protocol_id = find(ismember(tokens, "brainstorm_db"), 1) + 1;
protocol_name = tokens(protocol_id);

% Start brainstorm
if ~brainstorm('status')
    brainstorm nogui
end

% Open protocol
protocol = bst_get('Protocol', protocol_name);
gui_brainstorm('SetCurrentProtocol', protocol);

% Ensure FieldTrip is loaded
if ~exist('ft_defaults', 'file')
    [isInstalled, errMsg] = bst_plugin('Install', 'fieldtrip');
    if ~isInstalled
        error('Could not load FieldTrip plugin via Brainstorm: %s', errMsg);
    end
    bst_plugin('Load', 'fieldtrip');
end
ft_defaults;

%% ------------------------------------------------------------------
%% Step 1: Load brainstorm data and export to fieldtrip
%% ------------------------------------------------------------------

sFile = in_bst_data(DataFile);          % load file descriptor / raw link
ChannelFile= bst_get('ChannelFileForStudy', DataFile);
ChannelMat = in_bst_channel(ChannelFile);

if isstruct(sFile.F)
 % If raw is not brainstorm binary, stop
    if ~strcmp(sFile.F.format, "BST-BIN")
        error("Raw data files are only supported in brainstorm binary format (BST-BIN).");
    end

    % Load raw source
    [~, time, ~] = in_fread(sFile.F, channels, 1, []);
    
    % Grab BAD events for later
    bad_event_idx = find(contains({sFile.F.events.label}, 'BAD'));
 
else
    time = sFile.Time;

    % Grab BAD events for later
    bad_event_idx = find(contains({sFile.Events.label}, 'BAD'));

end

fs = round(1/(time(2) - time(1)));

% Remove bad meg channels/Only keep good channels
extracranialIdx = extracranialIdx(sFile.ChannelFlag(extracranialIdx) == 1);

if data_cov 
    allChanIdx = [extracranialIdx(:); intracranialIdx];
elseif noise_cov
    allChanIdx = [extracranialIdx(:)];
end

ftData = out_fieldtrip_data(DataFile, ChannelMat, allChanIdx, 0);
ftData.time = {time};

% ftData now has fields: .trial{1} [nChan x nSamples], .time{1}, .label, .fsample
fprintf('Exported %d channels, %d samples (%.1f min) at %d Hz\n', ...
    numel(ftData.label), numel(ftData.time{1}), ...
    numel(ftData.time{1})/fs/60, fs);

% remove bad segments

keepSamples = true(size(ftData.time{1}));
for iBad = 1:numel(bad_event_idx)
    for iEvents = 1:numel(sFile.F.events(bad_event_idx(iBad)).epochs)
        onset_bad = sFile.F.events(bad_event_idx(iBad)).times(1, iEvents);
        offset_bad = sFile.F.events(bad_event_idx(iBad)).times(2, iEvents);
        bad_idx = ftData.time{1} >= onset_bad & ftData.time{1} <= offset_bad;

        keepSamples(bad_idx) = false;
    end
end

% Baselinecorrect data, to not introduce concatination artifacts
CC = bwconncomp(keepSamples);
data_clean = zeros(size(ftData.trial{1}, 1), sum(keepSamples));
start_idx = 1;
for i=1:length(CC.PixelIdxList)
    end_idx = length(CC.PixelIdxList{i});
    data_clean(:, start_idx:start_idx + end_idx -1) = ftData.trial{1}(:, CC.PixelIdxList{i}) - mean(ftData.trial{1}(:, CC.PixelIdxList{i}), 2);
    start_idx = start_idx + end_idx;
end

% Update ftData Object
ftData.trial = {data_clean};
% Updating time to be continuous, otherwise the sampling rate is estimated
% wrong after epoching (step 2) -> start from 0s
ftData.time = {(0:(1/fs)*1000:(size(data_clean, 2)/fs * 1000)-1)/1000};

fprintf('Removed %d samples (%.1f min) due to marked bad segments.\n', ...
    sum(~keepSamples), ...
    sum(~keepSamples)/fs/60);

%% ------------------------------------------------------------------
%% Step 2: Cut continuous data into pseudo-trials for CSD estimation
%% ------------------------------------------------------------------

cfg = [];
cfg.length  = seg_length_sec;
cfg.overlap = seg_overlap;
ftEpoched = ft_redefinetrial(cfg, ftData);

n_segments = numel(ftEpoched.trial);
fprintf('Segmented into %d pseudo-trials of %.1f s each (%.0f%% overlap)\n', ...
    n_segments, seg_length_sec, seg_overlap*100);

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
cfg.channel    = ftEpoched.label;          % all MEG + LFP channels
cfg.channelcmb = {ftEpoched.label, ftEpoched.label};  % all pairwise combinations

freq_csd = ft_freqanalysis(cfg, ftEpoched);

fprintf('Computed CSD for band %s (center %.1f Hz, smoothing +/-%.1f Hz)\n', ...
    freq_id, centerFreq, halfBW);

% Save multilayer results
ResultsMat = db_template('noisecovmat');
ResultsMat.NoiseCov = freq_csd;

if data_cov
    OutputFile = fullfile(fileparts(DataFile), ...
        sprintf('ndatacov_CSD_fieldtrip_%s.mat', freq_id));
    covariance_comment = sprintf('Data Cross Spectral Density (%s, Fieldtrips cfg)', freq_id);
elseif noise_cov
    OutputFile = fullfile(fileparts(DataFile), ...
        sprintf('noisecov_CSD_fieldtrip_%s.mat', freq_id));
    covariance_comment = sprintf('Noise Cross Spectral Density (%s, Fieldtrips cfg)', freq_id);
end

ResultsMat.Comment = covariance_comment;


bst_save(OutputFile, ResultsMat, 'v6');

fprintf('Saved %s FieldTrip CSD to: %s\n', freq_id, OutputFile);

