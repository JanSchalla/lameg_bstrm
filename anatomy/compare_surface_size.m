function newTessFiles = compare_surface_size(TessFiles)

if ~iscell(TessFiles)
    error('Input must be of format "cell".');
end

% How many meshes/layers are combined
layers = numel(TessFiles);

% Sort files based on size -> Thus downsampling is minimized
size_info = cell(layers, 2);
for i = 1:layers
    VarInf = whos('-file', TessFiles{i}, 'Vertices');
    size_info{i, 1} = TessFiles{i};
    size_info{i, 2} = VarInf.size(1);
end

% Identify minimum vertices present on all surfaces
min_verts = min(size_info{:, 2});

newTessFiles = cell(layers, 1);
for i=1:layers
    dsTessFile = size_info{i, 1};
    if min_verts ~= size_info{i, 2}
        % Downsample bigger surfaces
        fprintf('Downsampling %s to %i Vertices.\n', size_info{i, 1}, min_verts);

        [dsTessFile, ~, ~, ~] = tess_downsize(size_info{i, 1}, min_verts, 'reducepatch');
        dsTessFile = file_fullpath(dsTessFile);
    end
    newTessFiles{i} = dsTessFile;
end