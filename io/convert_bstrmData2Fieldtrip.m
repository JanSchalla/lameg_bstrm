function ftData = convert_bstrmData2Fieldtrip(DataFile)

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
    [~, time, ~] = in_fread(sFile.F, ChannelMat, 1, []);
    
    % Grab BAD events for later
    bad_event_idx = find(contains({sFile.F.events.label}, 'BAD'));
 
else
    time = sFile.Time;

    % Grab BAD events for later
    bad_event_idx = find(contains({sFile.Events.label}, 'BAD'));

end

fs = round(1/(time(2) - time(1)));

% Remove bad meg channels/Only keep good channels
good_channel = find(sFile.ChannelFlag == 1);

ftData = out_fieldtrip_data(DataFile, ChannelMat, good_channel, 0);
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