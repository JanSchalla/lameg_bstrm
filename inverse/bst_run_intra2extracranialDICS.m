function bst_run_intra2extracranialDICS(DataFile, HeadModelFile, extracranialIdx, intracranialIdx, params)

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

% Get anatomy (Multilayer and one surface where the vertices are
% corresponding to the other).
multilayer_idx = find(contains({sSubject.Surface.Comment}, 'multilayer'));
corresponding_surf_idx = find(contains({sSubject.Surface.Comment}, 'corresponding'));

sFile = in_bst_data(DataFile);          % load file descriptor / raw link
ChannelFile= bst_get('ChannelFileForStudy', DataFile);
ChannelMat = in_bst_channel(ChannelFile);

fs = round(sFile.F.header.sfreq);

% Remove bad meg channels/Only keep good channels
extracranialIdx = extracranialIdx(sFile.ChannelFlag(extracranialIdx) == 1);
allChanIdx = [extracranialIdx(:); intracranialIdx];

ftData = out_fieldtrip_data(DataFile, ChannelMat, allChanIdx, 0);
ftData.time = {(sFile.Time)};

% ftData now has fields: .trial{1} [nChan x nSamples], .time{1}, .label, .fsample
fprintf('Exported %d channels, %d samples (%.1f min) at %d Hz\n', ...
    numel(ftData.label), numel(ftData.time{1}), ...
    numel(ftData.time{1})/fs/60, fs);

% remove bad segments
bad_event_idx = find(contains({sFile.F.events.label}, 'BAD'));
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
end

%% ------------------------------------------------------------------
%% Step 4: Load precomputed leadfield/head model from Brainstorm
%% ------------------------------------------------------------------

[ftHeadmodel, ftSourcemodel] = out_fieldtrip_headmodel(HeadModelFile, ChannelMat, extracranialIdx, 1);

nSources = numel(ftSourcemodel.leadfield);
nSurface = nSources / 2;
fprintf('Loaded leadfield for %d source points (%d per layer)\n', nSources, nSurface);

%% ------------------------------------------------------------------
%% Step 5: Run DICS with LFP as reference channel, per band
%% ------------------------------------------------------------------

dics_result = struct();

for b = 1:numel(bandNames)
    bname = bandNames{b};
    band = bandDefs.(bname);
    centerFreq = mean(band);

    cfg = [];
    cfg.method     = 'dics';
    cfg.refchan    = refchanLabel;
    %cfg.channel    = ftEpoched.label(~ismember(ftEpoched.label, ChannelMat.Channel(iLFP).Name));
    %cfg.frequency  = centerFreq;
    cfg.sourcemodel = ftSourcemodel;
    cfg.headmodel   = ftHeadmodel;        % may be unused if leadfield is precomputed, but FieldTrip still expects the field
    cfg.dics.projectnoise = 'yes';        % ???
    cfg.dics.lambda        = '5%';        % regularization, matches Brainstorm/FieldTrip common default
    %cfg.dics.keepfilter    = 'yes';       % needed if you want spatial filters for later time-domain projection
    %cfg.dics.realfilter    = 'yes';       % ???
    cfg.reducerank = 2;
    cfg.grad = freq_csd.(bname).grad;

    dics_result.(bname) = ft_sourceanalysis(cfg, freq_csd.(bname));

    fprintf('DICS complete for band %s\n', bname);
end

%% ------------------------------------------------------------------
%% Step 6: Extract coherence with the LFP reference, split by layer
%% ------------------------------------------------------------------

coh_white = struct();
coh_pial  = struct();
contrast  = struct();
layer_diff = struct();
epsilon = 1e-6;

for b = 1:numel(bandNames)
    bname = bandNames{b};
    coh_all = dics_result.(bname).avg.coh(:);   % [nSources x 1]

    coh_white.(bname) = coh_all(1:nSurface);
    coh_pial.(bname)  = coh_all(nSurface+1:end);

    layer_diff.(bname) = coh_pial.(bname) - coh_white.(bname);

    contrast.(bname) = log((coh_pial.(bname) + epsilon) ./ (coh_white.(bname) + epsilon));

    fprintf('Band %s: mean coh white=%.4f, pial=%.4f\n', bname, ...
        mean(coh_white.(bname)), mean(coh_pial.(bname)));
end


%% ------------------------------------------------------------------
%% Step 7: Export results back to Brainstorm for visualization
%% ------------------------------------------------------------------

for b = 1:numel(bandNames)
    bname = bandNames{b};

    % Save multilayer results
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = [coh_white.(bname); coh_pial.(bname)];  % full multilayer vector
    ResultsMat.Time          = 0;
    ResultsMat.Comment       = sprintf('DICS_coh_STN_%s', bname);
    ResultsMat.nComponents   = 1;
    ResultsMat.SurfaceFile   = sSubject.Surface(multilayer_idx).FileName;

    OutputFile = fullfile(fileparts(HeadModelFile), ...
        sprintf('results_DICS_coh_STN_%s.mat', bname));
    bst_save(OutputFile, ResultsMat, 'v6');

    fprintf('Saved %s DICS coherence map to: %s\n', bname, OutputFile);

    % Save layer fraction (logarythmic scale)
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = contrast.(bname);
    ResultsMat.Time          = 0;
    ResultsMat.Comment       = sprintf('DICS_coh_STN_fract_%s', bname);
    ResultsMat.nComponents   = 1;
    ResultsMat.SurfaceFile   = sSubject.Surface(corresponding_surf_idx(1)).FileName;

    OutputFile = fullfile(fileparts(HeadModelFile), ...
        sprintf('results_DICS_coh_STN_fract_%s.mat', bname));
    bst_save(OutputFile, ResultsMat, 'v6');

    fprintf('Saved %s DICS coherence fraction map to: %s\n', bname, OutputFile);

    % Save layer difference
    ResultsMat = db_template('resultsmat');
    ResultsMat.ImagingKernel = [];
    ResultsMat.ImageGridAmp  = layer_diff.(bname);
    ResultsMat.Time          = 0;
    ResultsMat.Comment       = sprintf('DICS_coh_STN_diff_%s', bname);
    ResultsMat.nComponents   = 1;
    ResultsMat.SurfaceFile   = sSubject.Surface(corresponding_surf_idx(1)).FileName;

    OutputFile = fullfile(fileparts(HeadModelFile), ...
        sprintf('results_DICS_coh_STN_diff_%s.mat', bname));
    bst_save(OutputFile, ResultsMat, 'v6');

    fprintf('Saved %s DICS coherence difference map to: %s\n', bname, OutputFile);
end
