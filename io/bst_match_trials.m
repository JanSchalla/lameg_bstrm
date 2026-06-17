function bst_match_trials(study_path, event_name1, event_name2)

study_mat = load(study_path);

nOrig_bad = numel(study_mat.BadTrials);

files = study_mat.BadTrials(:);

% Match event1 and event 2 to be bad if the other is bad
for i = 1:numel(files)
    if contains(files{i}, event_name1)
        opposite = strrep(files{i}, event_name1, event_name2);
    elseif contains(files{i}, event_name2)
        opposite = strrep(files{i}, event_name2, event_name1);
    else
        continue
    end

    if ~ismember(opposite, files)
        files{end+1} = opposite;
    end
end

nNew_bad = numel(files);

fprintf('Added %i Bad trials to match.', nNew_bad - nOrig_bad)

% Update brainstormstudy.mat
study_mat.BadTrials = files';

save(study_path, '-struct', 'study_mat');