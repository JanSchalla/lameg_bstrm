function bst_fT_calculateCSDMatrix(DataFile, extracranialIdx, intracranialIdx, params)

%% ------------------------------------------------------------------
%% Defaults
%% ------------------------------------------------------------------

seg_length_sec = 2;
seg_overlap = 0.5;
bandDefs = struct( ...
    'theta', [4, 7], ...
    'alpha', [7, 13], ...
    'beta', [13 30], ...
    'gamma', [30 60]);
multilayer = false;

if exist('params', 'var') && ~isempty(params)
    if isfield(params, 'seg_length_sec')
        seg_length_sec = params.seg_length_sec;
    end

    if isfield(params, 'seg_overlap')
        seg_overlap = params.seg_overlap;
    end

    if isfield(params, 'bandDefs')
        bandDefs = params.bandDefs;
    end

    if isfield(params, 'multilayer')
        multilayer = params.multilayer;
    end
end

% Verify input
bandNames = fieldnames(bandDefs);
for i = 1:numel(bandNames)
    if numel(bandDefs.(bandNames{i})) ~= 2
        error('Frequency band definiton needs one lower and one upper bound!');
    end
end

if seg_overlap < 0 || seg_overlap > 0.99
    error('Window overlap is not possible. Specify between 0 - 0.99.');
end

if multilayer
    warning('Multilayer analysis is set to TRUE.');
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

subject_id = find(ismember(tokens, "data"), 1) + 1;
subject_name = convertStringsToChars(tokens(subject_id));
study_name = convertStringsToChars(tokens(subject_id + 1));
data_name = convertStringsToChars(tokens(subject_id + 2));

% Start brainstorm
if ~brainstorm('status')
    brainstorm nogui
end

% Open protocol
protocol = bst_get('Protocol', protocol_name);
gui_brainstorm('SetCurrentProtocol', protocol);

sProtocol = bst_get('ProtocolSubjects');

subj_idx = find(ismember({sProtocol.Subject.Name}, subject_name));

sSubject = sProtocol.Subject(subj_idx);
clear("sProtocol");

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
allChanIdx = [extracranialIdx(:); intracranialIdx];

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

refchanLabel = ChannelMat.Channel(intracranialIdx).Name;

freq_csd = struct();

for b = 1:numel(bandNames)
    bname = bandNames{b};
    band = bandDefs.(bname);
    centerFreq = mean(band);
    halfBW = diff(band)/2;

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

    freq_csd.(bname) = ft_freqanalysis(cfg, ftEpoched);

    fprintf('Computed CSD for band %s (center %.1f Hz, smoothing +/-%.1f Hz)\n', ...
        bname, centerFreq, halfBW);

    % Reorder Cross-Spectrum into matrix (drop LFP
    chanlabels = freq_csd.(bname).label(~(contains(freq_csd.(bname).label, 'LFP')));     % canonical channel order
    Nchan = numel(chanlabels);
    data_cov = nan(Nchan, Nchan);

    [a_test, b_test] = build_csd_matrix(freq_csd.(bname),  ChannelMat.Channel(extracranialIdx).Name);
end
