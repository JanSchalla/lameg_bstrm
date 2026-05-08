function New_multilayer_fname = tess_create_multilayer_mesh(n_layers, wm_high_res, pial_high_res, params)

%% TESS_CREATE_MULTILAYER_MESH  Create Brainstorm-compatible multilayer cortical surface
%
%  DESCRIPTION:
%       Generates a multilayer cortical mesh with user-specified number of 
%       layers between white matter (wm_high_res) and pial (pial_high_res) 
%       surfaces. Intermediate surface are generated with linear spacing 
%       inbetween wm and pial surface, with optional non-linear spacing. 
%       Downsamples pial surfaces to matching vertex count, creates 
%       intermediate surfaces (both with corresponding vertices to the 
%       white matter surface), computes radial link vectors as vertex 
%       normals, concatenates layers, and updates scout atlases for 
%       Brainstorm visualization/processing compatibility.
%
%  INPUTS:
%       n_layers     - [scalar] Number of surface layers (including wm + pial)
%       wm_high_res  - [string] Full path to high-res white matter surface file (*.mat)
%       pial_high_res- [string] Full path to high-res pial surface file (*.mat) 
%       params       - [struct] (optional) Processing parameters:
%                       - .keep_proc_files  - [bool] Keep intermediate files (default: true)
%                       - .newNbVertices    - [scalar] Target vertex count (default: 15002)
%                       - .inflation_function- ['linear'|'cos'] Spacing function (default: 'linear')
%
%  OUTPUTS:
%       New_multilayer_fname - [string] Full path to created multilayer surface file
%
%  DEPENDENCIES:
%       - Brainstorm toolbox functions: tess_downsize(), in_tess_bst(), tess_concatenate(), 
%         bst_save(), bst_progress(), db_add_surface()
%
%  EXAMPLE:
%       params.linear = true;
%       params.newNbVertices = 15000;
%       fname = tess_create_multilayer_mesh(10, 'wm.mat', 'pial.mat', params);
%
%  AUTHORS:
%       [Jan Schalla], [Heinrich-Heine University], [2025]
%       Based on Brainstorm tessellation workflows
%
%  CHANGELOG:
%       [2025-09-26] - Initial implementation
%       [2026-02-16] - Expansion to create arbitrary amount of layers

%defaults
keep_proc_files = true;
newNbVertices = 15002;
inflation_function = 'linear';

% parse defaults
if exist('params', 'var') && ~isempty(params)
        if isfield(params, 'keep_proc_files')
            keep_proc_files = params.keep_proc_files;
        end

        if isfield(params, 'newNbVertices')
            newNbVertices = params.newNbVertices;
        end

        if isfield(params, 'inflation_function')
            inflation_function = params.inflation_function;
        end
end


%% User feedback
fprintf('A %i-layer mesh will be created.\n', n_layers);
fprintf('Inner boundary will be %s.\n', wm_high_res);
fprintf('Outer bounday will be %s.\n', pial_high_res);

if strcmp(inflation_function, 'linear')
    inflation_steps = linspace(0, 1, n_layers);
elseif strcmp(inflation_function, 'cos')
    inflation_steps = (1-cos(linspace(0,pi,n_layers)))/2;
end


if isnumeric(newNbVertices)
    if newNbVertices == floor(newNbVertices)
        fprintf('Nr. of vertices to downsample specified to %i\n', ...
            newNbVertices);
    else
        fprintf('Nr. of specified verticies not an integer. Rounding down to %i\n', ...
            floor(newNbVertices));
        newNbVertices = floor(newNbVertices);
    end
else
    error('Specified Nr. of vertices is not numeric.');
end

bst_progress('start', 'Create Multilayer mesh', 'Processing... ');

%% Downsample
% Downsample first surface
% I corresponds to the vertices/rows which are in the high-res surface and also
% in the downsampled surface
% J corresponds to the vertices, which are in the downsampled surface and
% also in the high-res surface

fprintf('Downsampling white matter surface to %i Vertices.\n', newNbVertices);
[ds_firstSurf, ~, I, ~] = tess_downsize(wm_high_res, newNbVertices, 'reducepatch');
vertices_per_surface = numel(I);

ds_firstSurf = file_fullpath(ds_firstSurf);

% Ask the user if scouts should be created
resp = questdlg('Do you want to create individual scouts now?', ...
                'Create Scouts', ...
                'Yes','No','No');

if strcmp(resp,'Yes')
    % Open the surface in Brainstorm viewer
    hFig = view_surface(ds_firstSurf);
    bst_figures('SetCurrentFigure', hFig);     % optional, make it active
    
    % Optional message to guide the user
    disp('Surface opened. Use Brainstorm GUI: Scouts > Create scout.');
        
    % Ask the user if scouts should be created
    disp('Surface opened. Create/edit scouts in Brainstorm, then press any key here to continue...');
    pause;    % waits for key press in command window
    disp('Continuing script after scout editing (Make sure to close the brainstorm figure to copy scouts to other layers!.');
end

% Load in downsampled surface
TessMat_wm_ds = in_tess_bst(ds_firstSurf);

% % Update naming & comment of the first surface
% tok = regexp(ds_firstSurf, '(\d+)V\.mat', 'tokens');
% nVerts_path = str2num(tok{1}{1});

TessMat_wm_ds.Comment = sprintf('%s_corresponding', TessMat_wm_ds.Comment);
% nVerts_new = size(TessMat_wm_ds.Vertices, 1);

% ds_firstSurf_new = strrep(ds_firstSurf, sprintf('%dV.mat', nVerts_path), sprintf('%dV.mat', nVerts_new));

bst_save(ds_firstSurf, TessMat_wm_ds, 'v7');
% Delete old file
% delete(ds_firstSurf);
% 
% % Overwrite path
% ds_firstSurf = ds_firstSurf_new;

layer_fnames = cell(n_layers, 1);

layer_fnames{1} = ds_firstSurf;

%% Resample rest of surfaces with the information from first surface
% First resample pial surface
TessMat_pial_ds = in_tess_bst(pial_high_res);
TessMat_pial_ds.Vertices = TessMat_pial_ds.Vertices(I, :);
%Copy Structures
TessMat_pial_ds.Faces = TessMat_wm_ds.Faces;
TessMat_pial_ds.VertConn = TessMat_wm_ds.VertConn;
TessMat_pial_ds.Curvature = TessMat_wm_ds.Curvature;
TessMat_pial_ds.SulciMap = TessMat_wm_ds.SulciMap;
TessMat_pial_ds.Atlas = TessMat_wm_ds.Atlas;
TessMat_pial_ds.VertNormals = TessMat_wm_ds.VertNormals;
TessMat_pial_ds.Reg = TessMat_wm_ds.Reg;
% Save
[f_path,f_name, ~] = fileparts(file_fullpath(pial_high_res));
tokens = split(f_name, '_');
surf_ident_idx = find(ismember(tokens, 'cortex')) + 1;
%Update comments
TessMat_pial_ds.Comment = sprintf('%s_%iV_corresponding', tokens{surf_ident_idx}, size(TessMat_pial_ds.Vertices, 1));
dest_ds_Surf = fullfile(f_path, sprintf('tess_cortex_%s_%iV_ds_correspondingVerts.mat', tokens{surf_ident_idx}, length(I)));
idx = 1;
while isfile(dest_ds_Surf)
    fprintf('Renaming %s, due to already being present in database.\n', dest_ds_Surf);
    idx = idx + 1;
    dest_ds_Surf = fullfile(f_path, sprintf('tess_cortex_%s_%iV_ds_correspondingVerts_0%i.mat', tokens{surf_ident_idx}, length(I), idx));
end
    
bst_save(dest_ds_Surf, TessMat_pial_ds, 'v7');
layer_fnames{end} = dest_ds_Surf;

for s=2:n_layers-1
    inflation_perc = inflation_steps(s);
    fprintf('Processing surface %i of %i. Creating intermediate surface (%.3f)...\n', s, n_layers-1,inflation_perc); 
    TessMat_ds = TessMat_wm_ds;
        
    % Take vertices present in first surface to keep in second surface
    TessMat_ds.Vertices =  TessMat_wm_ds.Vertices + inflation_perc * (TessMat_pial_ds.Vertices - TessMat_wm_ds.Vertices);
    % Copy structs
    TessMat_ds.Faces = TessMat_wm_ds.Faces;
    TessMat_ds.VertConn = TessMat_wm_ds.VertConn;
    TessMat_ds.Curvature = TessMat_wm_ds.Curvature;
    TessMat_ds.SulciMap = TessMat_wm_ds.SulciMap;
    TessMat_ds.Atlas = TessMat_wm_ds.Atlas;
    TessMat_ds.VertNormals = TessMat_wm_ds.VertNormals;
    TessMat_ds.Reg = TessMat_wm_ds.Reg;

    % Save
    % Update comments
    TessMat_ds.Comment = sprintf('%.3f_%iV', inflation_perc, size(TessMat_ds.Vertices, 1));
    
    dest_ds_Surf = fullfile(f_path, sprintf('tess_cortex_%.3f_%iV_ds_correspondingVerts.mat', inflation_perc, length(I)));
    idx = 1;
    while isfile(dest_ds_Surf)
        fprintf('Renaming %s, due to already being present in database.\n', dest_ds_Surf);
        idx = idx + 1;
        dest_ds_Surf = fullfile(f_path, sprintf('tess_cortex_%s_%iV_ds_correspondingVerts_0%i.mat', tokens{surf_ident_idx}, length(I), idx));
    end
        
    bst_save(dest_ds_Surf, TessMat_ds, 'v7');
    layer_fnames{s} = dest_ds_Surf;
end

%% Compute Normals 
% Compute link vectors 
% (Work best when doing laMEG: 10.1016/j.neuroimage.2020.116862)
% TO-DO: Offer up alternative ways of computing link vectors

fprintf('Computing link vectors between smallest and biggest surface ...\n');
small_verts = TessMat_wm_ds.Vertices;
big_verts = TessMat_pial_ds.Vertices;

if size(small_verts, 1) ~= size(big_verts, 1)
    error('Surfaces must have the same number of vertices. Check Downsampling. Aborting ...');
end

link_vectors = small_verts - big_verts;
norm_e = sqrt(sum(abs(link_vectors).^2, 2)); %Normalize by euclidian norm
norm_e(norm_e < eps) = 1; % Set to small norms to 1 to avoid computaional problems later on.

link_vectors = link_vectors ./ norm_e; % Bring link vectors to unit length.

%% Update VertNormals to be link vectors
for iFile = 1:length(layer_fnames)
    TessMat_ds = in_tess_bst(layer_fnames{iFile});
    fprintf('Updating link vectors for %s\n', layer_fnames{iFile});
    if all(size(TessMat_ds.VertNormals) == size(link_vectors))
        TessMat_ds.VertNormals = link_vectors;
    elseif all(size(TessMat_ds.VertNormals) == size(link_vectors'))
        TessMat_ds.VertNormals = link_vectors';
        disp('Link vectors have been transposed to match size of original .VertNormals.');
    else
        error('Size of original VertNormals and link vectors do not match! Aborting ...');
    end
    bst_save(layer_fnames{iFile}, TessMat_ds, 'v7');
end

%% Combine Layers
fprintf('Creating Multilayer surface using tess_concatenate ...\n');
out_name = sprintf('tess_cortex_multilayer_%i', n_layers);

% Concatate the n-layers to one surface
multilayer_fname = tess_concatenate(layer_fnames, out_name);
multilayer_fname = file_fullpath(multilayer_fname);

[anat_path, ~, multi_ext] = fileparts(multilayer_fname);
New_multilayer_fname = sprintf('%s/%s%s', anat_path, out_name, multi_ext);
movefile(multilayer_fname, New_multilayer_fname);

%% Clean up brainstorm db from intermediate files
% Move created file so brainstorm does not get confused
if keep_proc_files

    disp('Renaming downsampled files.');

    for iFile = 1:length(layer_fnames)
        [anat_path, f_name, ~] = fileparts(layer_fnames{iFile}); 

        if iFile == 1
            tokens = split(f_name, '_');
            surf_ident_idx = find(ismember(tokens, 'cortex')) + 1;
            tmp_ds_surf = fullfile(anat_path, sprintf('tess_cortex_%s_%iV_ds_correspondingVerts.mat', tokens{surf_ident_idx}, length(I)));
        else
            tmp_ds_surf = fullfile(anat_path, sprintf('%s.mat', f_name));
        end

        if strcmp(layer_fnames{iFile}, tmp_ds_surf)
            continue
        else
            movefile(layer_fnames{iFile}, tmp_ds_surf)
        end
    end
    % Do nothing!
    % tmp_folder = fullfile(anat_path, 'tmp_multilayer');
    % fprintf('Moving downsampled files to %s ...\n', tmp_folder);
    % if not(isfolder(tmp_folder))
    %     mkdir(tmp_folder);
    % end
    % 
    % for iFile = 1:length(layer_fnames)
    %     [~, f_name, ~] = fileparts(layer_fnames{iFile}); 
    % 
    %     if iFile == 1
    %         tokens = split(f_name, '_');
    %         surf_ident_idx = find(ismember(tokens, 'cortex')) + 1;
    %         tmp_ds_surf = fullfile(tmp_folder, sprintf('tess_cortex_%s_%iV_ds_correspondingVerts.mat', tokens{surf_ident_idx}, length(I)));
    %     else
    %         tmp_ds_surf = fullfile(tmp_folder, sprintf('%s.mat', f_name));
    %     end
    %     movefile(layer_fnames{iFile}, tmp_ds_surf)
    % end
else
    fprintf('Deleting downsampled files ...\n');
    for iFile = 1:length(layer_fnames)
        delete(layer_fnames{iFile});
    end
end

%% Update scouts to be multilayer
TessMat = in_tess_bst(New_multilayer_fname);
fprintf('Updating multilayer scouts ...\n');
new_atlas_struct = db_template('atlas');
for i = 1:length(TessMat.Atlas)
    unique_scouts = unique({TessMat.Atlas(i).Scouts.Label});
    new_atlas_struct(i).Scouts = repmat(db_template('scout'), 0);
    new_atlas_struct(i).Name = TessMat.Atlas(i).Name;
    for ii=1:length(unique_scouts)
        duplicate_scouts = find(ismember({TessMat.Atlas(i).Scouts.Label}, unique_scouts{ii}));
        new_atlas_struct(i).Scouts(ii).Vertices = [TessMat.Atlas(i).Scouts(duplicate_scouts).Vertices]; % Combine vertices from both surfaces together
        %Take values from first original scout
        new_atlas_struct(i).Scouts(ii).Seed = TessMat.Atlas(i).Scouts(duplicate_scouts(1)).Seed;
        new_atlas_struct(i).Scouts(ii).Color = TessMat.Atlas(i).Scouts(duplicate_scouts(1)).Color;
        new_atlas_struct(i).Scouts(ii).Label = TessMat.Atlas(i).Scouts(duplicate_scouts(1)).Label;
        new_atlas_struct(i).Scouts(ii).Function = TessMat.Atlas(i).Scouts(duplicate_scouts(1)).Function;
        new_atlas_struct(i).Scouts(ii).Region = TessMat.Atlas(i).Scouts(duplicate_scouts(1)).Region;
        new_atlas_struct(i).Scouts(ii).Handles = TessMat.Atlas(i).Scouts(duplicate_scouts(1)).Handles;
    end

    if strcmp(TessMat.Atlas(i).Name, 'Structures')
        offsetVertices = 0;
        for iii = 1:n_layers
            if iii == 1
                label = 'Cortex white';
            elseif iii > 1 && iii < n_layers
                label = sprintf('Cortex %.3f', inflation_steps(iii));
            else
                label = 'Cortex pial';
            end
            new_atlas_struct(i).Scouts(end+1) = db_template('scout');
            % Create Info for single surfaces
            new_atlas_struct(i).Scouts(end).Vertices = 1+offsetVertices:vertices_per_surface+offsetVertices;
            new_atlas_struct(i).Scouts(end).Seed = max(new_atlas_struct(i).Scouts(end).Vertices);
            if iii == 1
                new_atlas_struct(i).Scouts(end).Color = [1.0 0.75 0.80]; % White matter scout in pastel pink
            elseif iii == n_layers
                new_atlas_struct(i).Scouts(end).Color = [0.3 0.75 0.93]; % Pial scout in light blue
            else
                new_atlas_struct(i).Scouts(end).Color = [1.0 0.65 0.4]; % Intermediate Surfaces in orange
            end
            new_atlas_struct(i).Scouts(end).Label = char(label);
            new_atlas_struct(i).Scouts(end).Region = 'UU';            
            
            fprintf('Multilayer sub-scout created for: \n');
            fprintf('Label: %s\n', label);
            fprintf('Seed: %i\n',new_atlas_struct(i).Scouts(end).Seed);
            fprintf('Start Vertex = %i, End Vertex = %i\n',1+offsetVertices, vertices_per_surface+offsetVertices);

            offsetVertices = offsetVertices + vertices_per_surface;
        end
    end
end

TessMat.Atlas = new_atlas_struct;
%If adding this somehow database gets corrupted and cannot open multilayer
%mesh
%TessMat = bst_history('add', TessMat, 'merge multilayer scouts', 'Merge completed');

bst_save(New_multilayer_fname, TessMat, 'v7');

disp('###################');
disp('Finishing Up ...');

%% Finish Up
% Make output filename relative
New_multilayer_fname_short = file_short(New_multilayer_fname);
% Get subject
[~, iSubject] = bst_get('SurfaceFile', wm_high_res);
% Register this file in Brainstorm database
[~] = db_add_surface(iSubject, New_multilayer_fname_short, out_name, 'Cortex');

bst_progress('stop');

disp('Done!');
fprintf('Multilayer surface is saved here: %s\n ', New_multilayer_fname_short)
disp('###################');