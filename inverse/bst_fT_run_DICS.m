function dics_result = bst_fT_run_DICS(DataFile, HeadModelFile, csd_cfg, params)

%% ------------------------------------------------------------------
%% Defaults
%% ------------------------------------------------------------------

multilayer = false;
sr_id = 'None'; % Identifier for SR
save_results = true;
iWhite = [];

if exist('params', 'var') && ~isempty(params)
    if isfield(params, 'multilayer')
        multilayer = params.multilayer;
    end

    if isfield(params, 'sr_id')
        sr_id = params.sr_id;
    end

    if isfield(params, 'save_results')
        save_results = params.save_results;
    end
    
    if isfield(params, 'iWhite')
        iWhite = params.iWhite;
    end
end

if strcmp(sr_id, 'None')
    warning('No Identifier for the source reconstruction specified. Saved results will not be identifiable by name.');
end

if ~save_results
    warning('Results are not saved!')
end

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

corresponding_surf_idx = find(contains({sSubject.Surface.Comment}, 'corresponding'));

%% ------------------------------------------------------------------
%% Step 1: Load brainstorm data and export to fieldtrip
%% ------------------------------------------------------------------

ChannelFile= bst_get('ChannelFileForStudy', DataFile);
ChannelMat = in_bst_channel(ChannelFile);

MEG_idx = find(contains({ChannelMat.Channel.Type}, 'MEG'));

% Update MEG_idx to only contain channels also present in cfg_csd!
good_chans = ismember({ChannelMat.Channel(MEG_idx).Name}, csd_cfg.label);
MEG_idx = MEG_idx(good_chans);

% % If not whitener supplied, multiply leadfield by identity -> same as doing
% % nothing
if isempty(iWhite)
    iWhite = eye(numel(MEG_idx));
end

%% Assumption: Supplied cfg_csd only contains MEG & one reference channel!
refchanLabel = csd_cfg.label{~contains(csd_cfg.label, 'MEG')};

%% ------------------------------------------------------------------
%% Step 4: Load precomputed leadfield/head model from Brainstorm
%% ------------------------------------------------------------------
bst_headmodel = in_bst_headmodel(HeadModelFile);

% Multiply whitener by leadfield to bring the leadfield into the same space
% as the provided crsspectrum -> if no whitener is specified (assuming no
% whitener applied), leadfield is multiplied by identity, leaving it as is.
bst_headmodel.Gain(MEG_idx, :) = iWhite * bst_headmodel.Gain(MEG_idx, :);

[ftHeadmodel, ftSourcemodel] = out_fieldtrip_headmodel(bst_headmodel, ChannelMat, MEG_idx, 1);

nSources = numel(ftSourcemodel.leadfield);
nSurface = nSources / 2;
fprintf('Loaded leadfield for %d source points (%d per layer)\n', nSources, nSurface);

%% ------------------------------------------------------------------
%% Step 5: Run DICS with LFP as reference channel, per band
%% ------------------------------------------------------------------

cfg = [];
cfg.method     = 'dics';
cfg.refchan    = refchanLabel;
%cfg.channel    = ftEpoched.label(~ismember(ftEpoched.label, ChannelMat.Channel(iLFP).Name));
%cfg.frequency  = centerFreq;
cfg.sourcemodel = ftSourcemodel;
cfg.headmodel   = ftHeadmodel;        % may be unused if leadfield is precomputed, but FieldTrip still expects the field
cfg.dics.projectnoise = 'yes';        % ???
cfg.dics.lambda        = '5%';        % regularization, matches Brainstorm/FieldTrip common default
cfg.dics.keepfilter    = 'yes';       % needed if you want spatial filters for later time-domain projection
cfg.dics.realfilter    = 'yes';       % ???
cfg.dics.fixedori      = 'yes';       % Needed to calcualte single orientation filter
cfg.dics.weightnorm    = 'nai';       % Needed to calcualte MNPSP
cfg.dics.keepcsd       = 'yes';       % Not sure if needed
cfg.reducerank         = 2;           % I get the warning, not used
cfg.grad = csd_cfg.grad;
cfg.whiten = 'yes';

dics_result = ft_sourceanalysis(cfg, csd_cfg);

fprintf('DICS complete for band %s\n', sr_id);

%% ------------------------------------------------------------------
%% Step 6: Extract coherence with the LFP reference, split by layer
%% ------------------------------------------------------------------

if multilayer

    epsilon = 1e-6;
    
    coh_all = dics_result.avg.coh(:);   % [nSources x 1]

    coh_white = coh_all(1:nSurface);
    coh_pial  = coh_all(nSurface+1:end);

    layer_diff = coh_pial - coh_white;

    contrast = log((coh_pial + epsilon) ./ (coh_white + epsilon));

    [wCoh, wCoh_idx] = max(coh_white);
    [pCoh, pCoh_idx] = max(coh_pial);

    fprintf('Band %s: max coh white=%.4f (Vertex: %i), pial=%.4f (Vertex: %i)\n', sr_id, ...
        wCoh, wCoh_idx, pCoh, pCoh_idx+nSurface);
else  
    coh_all = dics_result.avg.coh(:);
    
    fprintf('Band %s: mean coh=%.4f\n', sr_id, ...
        max(coh_all));
end

%% ------------------------------------------------------------------
%% Step 7: Export results back to Brainstorm for visualization
%% ------------------------------------------------------------------
if save_results
    if multilayer
        % Save multilayer results
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = [];
        ResultsMat.ImageGridAmp  = [coh_white; coh_pial];  % full multilayer vector
        ResultsMat.Time          = 0;
        ResultsMat.Comment       = sprintf('DICS_coh_%s_%s', refchanLabel, sr_id);
        ResultsMat.nComponents   = 1;
        ResultsMat.SurfaceFile   = bst_headmodel.SurfaceFile;
    
        OutputFile = fullfile(fileparts(HeadModelFile), ...
            sprintf('results_DICS_coh_%s_%s.mat', refchanLabel, sr_id));
        bst_save(OutputFile, ResultsMat, 'v6');
    
        fprintf('Saved %s DICS coherence map to: %s\n', sr_id, OutputFile);
    
        % % Save layer fraction (logarythmic scale)
        % ResultsMat = db_template('resultsmat');
        % ResultsMat.ImagingKernel = [];
        % ResultsMat.ImageGridAmp  = contrast;
        % ResultsMat.Time          = 0;
        % ResultsMat.Comment       = sprintf('DICS_coh_STN_fract_%s', sr_id);
        % ResultsMat.nComponents   = 1;
        % ResultsMat.SurfaceFile   = sSubject.Surface(corresponding_surf_idx(1)).FileName;
        % 
        % OutputFile = fullfile(fileparts(HeadModelFile), ...
        %     sprintf('results_DICS_coh_STN_fract_%s.mat', sr_id));
        % bst_save(OutputFile, ResultsMat, 'v6');
        % 
        % fprintf('Saved %s DICS coherence fraction map to: %s\n', sr_id, OutputFile);
        % 
        % % Save layer difference
        % ResultsMat = db_template('resultsmat');
        % ResultsMat.ImagingKernel = [];
        % ResultsMat.ImageGridAmp  = layer_diff;
        % ResultsMat.Time          = 0;
        % ResultsMat.Comment       = sprintf('DICS_coh_STN_diff_%s', sr_id);
        % ResultsMat.nComponents   = 1;
        % ResultsMat.SurfaceFile   = sSubject.Surface(corresponding_surf_idx(1)).FileName;
        % 
        % OutputFile = fullfile(fileparts(HeadModelFile), ...
        %     sprintf('results_DICS_coh_STN_diff_%s.mat', sr_id));
        % bst_save(OutputFile, ResultsMat, 'v6');
        % 
        % fprintf('Saved %s DICS coherence difference map to: %s\n', sr_id, OutputFile);
    
        % Save filter
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = vertcat(dics_result.avg.filter{:});
        ResultsMat.ImageGridAmp  = [];
        ResultsMat.Comment       = sprintf('DICS_coh_%s_FilterWeighhts_%s', refchanLabel, sr_id);
        ResultsMat.nComponents   = 1;
        ResultsMat.Time          = 0;
        ResultsMat.SurfaceFile   = bst_headmodel.SurfaceFile;
        
        OutputFile = fullfile(fileparts(HeadModelFile), ...
            sprintf('results_DICS_filter_%s_%s.mat', refchanLabel, sr_id));
        bst_save(OutputFile, ResultsMat, 'v6');
    
        fprintf('Saved %s DICS Filter to: %s\n', sr_id, OutputFile);
    else
        % Save multilayer results
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = [];
        ResultsMat.ImageGridAmp  = [coh_all];  % full multilayer vector
        ResultsMat.Time          = 0;
        ResultsMat.Comment       = sprintf('DICS_coh_STN_%s', sr_id);
        ResultsMat.nComponents   = 1;
        ResultsMat.SurfaceFile   = bst_headmodel.SurfaceFile;
    
        OutputFile = fullfile(fileparts(HeadModelFile), ...
            sprintf('results_DICS_coh_STN_%s.mat', sr_id));
        bst_save(OutputFile, ResultsMat, 'v6');
    
        fprintf('Saved %s DICS coherence map to: %s\n', sr_id, OutputFile);
    end
end