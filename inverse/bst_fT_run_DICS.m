function bst_fT_run_DICS(DataFile, HeadModelFile, csd_cfg, params)

%% ------------------------------------------------------------------
%% Defaults
%% ------------------------------------------------------------------

multilayer = false;

if exist('params', 'var') && ~isempty(params)
    if isfield(params, 'multilayer')
        multilayer = params.multilayer;
    end
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

ChannelFile= bst_get('ChannelFileForStudy', DataFile);
ChannelMat = in_bst_channel(ChannelFile);




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
    %band = bandDefs.(bname);
    %centerFreq = mean(band);

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
    cfg.reducerank         = 2;
    cfg.grad = freq_csd.(bname).grad;

    dics_result.(bname) = ft_sourceanalysis(cfg, freq_csd.(bname));

    fprintf('DICS complete for band %s\n', bname);
end

%% ------------------------------------------------------------------
%% Step 6: Extract coherence with the LFP reference, split by layer
%% ------------------------------------------------------------------

if multilayer

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
else
    for b = 1:numel(bandNames)
        bname = bandNames{b};
        
        coh_all.(bname) = dics_result.(bname).avg.coh(:);
        
        fprintf('Band %s: mean coh=%.4f\n', bname, ...
            mean(coh_all.(bname)));
    end
end

%% ------------------------------------------------------------------
%% Step 7: Export results back to Brainstorm for visualization
%% ------------------------------------------------------------------

if multilayer
    for b = 1:numel(bandNames)
        bname = bandNames{b};
        % Save multilayer results
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = dics_result.(bname).avg.filter;
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
    
        % Save filter and C matrices
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = vertcat(dics_result.(bname).avg.filter{:});
        ResultsMat.CSD_data = freq_csd.(bname) 
    end
end