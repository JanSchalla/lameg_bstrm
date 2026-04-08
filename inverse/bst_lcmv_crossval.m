function bst_lcmv_crossval(sFiles, nfolds)

crossval_results = struct();

bstrm_path = file_fullpath(sFiles{1});

if length(split(bstrm_path, "/")) > length(split(bstrm_path, "\"))
        separator = "/";
    else
        separator = "\";
end

tokens = split(bstrm_path, separator);
    
protocol_id = find(ismember(tokens, "brainstorm_db"), 1) + 1;
protocol_name = tokens(protocol_id);

subject_id = find(ismember(tokens, "data"), 1) + 1;
subject_name = convertStringsToChars(tokens{subject_id});
study_name = convertStringsToChars(tokens{subject_id + 1});
data_name = convertStringsToChars(tokens{subject_id + 2});

% Start brainstorm
if ~brainstorm('status')
    brainstorm nogui
end

% Open protocol
protocol = bst_get('Protocol', protocol_name);
gui_brainstorm('SetCurrentProtocol', protocol);

% Get relevant study
study_filename = sprintf('%s/%s/brainstormstudy.mat', subject_name, study_name);
[study, study_id] = bst_get('Study', study_filename);

% Get channels
channels = in_bst_channel(study.Channel.FileName);

data_filename = convertStringsToChars(sprintf("%s/%s/%s", subject_name, study_name, data_name));
data_struct = in_bst_data(data_filename);

% Get only MEG channels
iMeg = channel_find(channels.Channel, 'MEG');
% Get valid channels, not marked as bad prior to cross validation
valid_chans = iMeg(data_struct.ChannelFlag(iMeg) == 1);

%(save)
orig_bad_chans = iMeg(data_struct.ChannelFlag(iMeg) < 0);

% Total number of channels usable for SR
Nchans = length(valid_chans);

% Partition data into n-folds
fold_size = floor(Nchans/nfolds);

folds = cell(nfolds, 1);
pre_fold_chans = valid_chans;
for f=1:nfolds
    % % Iterative approach 
    % start_idx = (f-1)*fold_size + 1;
    % end_idx = min(f*fold_size, Nchans);
    % folds{f} = valid_chans(start_idx:end_idx);

    % Random approach
    % Take random channels into each fold
    idx = zeros(1, length(pre_fold_chans)) ~= 0;
    idx(1:fold_size) = true;
    
    idx = (idx(randperm(length(idx))));
    folds{f} = pre_fold_chans(idx);
    pre_fold_chans = pre_fold_chans(~idx);
end

% % Just a check if the random approach worked
% duplicate_chans = 0;
% for i = 1:nfolds
%     for j=1:nfolds
%         if j == i
%             continue
%         end
%         duplicate_chans = duplicate_chans + sum(folds{i} == folds{j});
%     end
% end
% disp(duplicate_chans);

% Pre-allocate output measures
crosserr = zeros(nfolds, fold_size);  % RMS error per test channel
allrms  = zeros(nfolds, fold_size);  % RMS signal per test channel
allfract = zeros(nfolds, fold_size); % Fractional RMS error

fprintf('Running %i-fold cross validation ...\n', nfolds);
for b=1:nfolds
    test_chans = folds{b};
    Ntest = length(test_chans);
    fprintf('Fold %i/%i\n', b, nfolds);

    %% Step 1: Mark test channels as bad
    % As source reconstruction in brainstorm works if the channels are bad in one trial they
    % are treated as bad for all trials -> good for us less computationally
    % heavy!
    data_struct = in_bst_data(sFiles{1});

    % Set test channels for the fold to bad
    fprintf('Setting channels as bad: \n%s, %s, %s, %s, %s, %s, %s, %s, %s, %s,\n%s, %s, %s, %s, %s, %s, %s, %s, %s, %s,\n%s, %s, %s, %s, %s, %s, %s, %s, %s\n', channels.Channel(test_chans(:)).Name)%, channels.Channel(folds{b}(2)).Name)
    data_struct.ChannelFlag(test_chans) = -1;

    % here alway sFiles{1} is altered to have the test channels marked as
    % bad.
    %% Be carefull when debugging to always change it back to the original
    %% valid_chans!!
    bst_save(file_fullpath(sFiles{1}), data_struct, 'v7', 0);

    %% Setp 2: Calculate LCMV Kernel only with training data
    results_file = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
        'output', 1, ...
        'inverse', struct(...
            'Comment',        char(sprintf('PNAI: MEG ALL (fold %i)', b)), ...
            'InverseMethod',  'lcmv', ...
            'InverseMeasure', 'nai', ...
            'SourceOrient',   {{'fixed'}}, ...
            'Loose',          0.2, ...
            'UseDepth',       0, ...
            'WeightExp',      0.5, ...
            'WeightLimit',    10, ...
            'NoiseMethod',    'median', ...
            'NoiseReg',       0.1, ...
            'SnrMethod',      'rms', ...
            'SnrRms',         1e-06, ...
            'SnrFixed',       3, ...
            'ComputeKernel',  1, ...
            'DataTypes',      {{'MEG GRAD', 'MEG MAG'}})); 

    %% Step 3: Restore originial bad channels
    % Set all to good
    data_struct.ChannelFlag(iMeg) = 1;
    % Apply original bad channels
    data_struct.ChannelFlag(orig_bad_chans) = -1; 
    bst_save(file_fullpath(sFiles{1}), data_struct, 'v7', 0);

    % Load imaging kernel (computed with training data)
    resultMat = in_bst_results(results_file(1).FileName, 0);

    %% Step 4: Get full leadfield (all MEG sensors)
    HeadModelMat = in_bst_headmodel(resultMat.HeadModelFile);
    HeadModelGain = HeadModelMat.Gain; %get gain for each source for each sensor in xyz orientation
    % Use this and Thomas way of bringing source data into 3 orientations
    % or use the Leadfield estimation from bstrms in process_inverse_2018?
    
    % Here the leadfield is calculated like brainstorm does it ->
    % nChannels*nSources
    Wq = cell(1,size(HeadModelMat.GridOrient, 1));
    for i = 1:size(HeadModelMat.GridOrient, 1)
        tmp = HeadModelMat.GridOrient(i, :)';
        Wq{i} = tmp/norm(tmp);
    end
    
    WQ = blkdiag(Wq{1},Wq{2:end});
    
    % unwhitend leadfield
    L = HeadModelGain * WQ;

    % which to use? white or non-whitened?
    % Whitening i can only do when i have the whitener calculated on all
    % channels once and use it here!
    % Lw = resultMat.DataWhitener * L;


    %% STEP 5: Per-trial prediction and error computation
    errpred = 0;
    rmstot  = 0;
    fracterr = 0;

    for trial = 1:length(sFiles)
        % load original sensor data 
        data_struct_orig = in_bst_data(sFiles{trial});

        %reconstruct sensor data based on training kernel
        source_training = resultMat.ImagingKernel * data_struct_orig.F(resultMat.GoodChannel, :);

        % Backproject to all sensors
        y_pred_full = L * source_training;
        % Take only test channels 
        y_pred = y_pred_full(test_chans, :);

        err2 = (y_pred - data_struct_orig.F(test_chans, :)).^2;
        sig2 = data_struct_orig.F(test_chans, :).^2;
        diff2 = err2./sig2;
        rmserr = sqrt(mean(err2, 2));
        rmssig = sqrt(mean(sig2, 2));
        rmsdiff = sqrt(mean(diff2, 2));

        errpred = errpred + rmserr;
        rmstot = rmstot + rmssig;
        fracterr = fracterr + rmsdiff;
    end
    crosserr(b, :) = errpred./length(sFiles);
    allrms(b, :) = rmstot./length(sFiles);
    allfract(b, :) = fracterr./length(sFiles);
end

crossval_results.crosserr = crosserr;
crossval_results.allrms = allrms;
crossval_results.allfract = allfract;
crossval_results.test_chans_per_fold = folds;



h = plot(crossval_results.crosserr) %Error between the predicted and measured data (squared)
rowNames = arrayfun(@(i) sprintf('Fold-%d', i), 1:10, 'UniformOutput', false);
legend(h, rowNames, 'Location', 'best');

h = plot(crossval_results.allrms) %squared amplitude of measured signal
rowNames = arrayfun(@(i) sprintf('Fold-%d', i), 1:10, 'UniformOutput', false);
legend(h, rowNames, 'Location', 'best');

h = plot(crossval_results.allfract) % difference between the actual measured amplitude and the error
rowNames = arrayfun(@(i) sprintf('Fold-%d', i), 1:10, 'UniformOutput', false);
legend(h, rowNames, 'Location', 'best');