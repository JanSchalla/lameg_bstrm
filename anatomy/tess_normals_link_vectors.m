function multilayer_struct = tess_normals_link_vectors(multilayer_path)
% This works currently only for 2-layer surfaces!

[multilayer_struct, ~] = in_tess_bst(multilayer_path);

small_surf_idx = 1:size(multilayer_struct.Vertices, 1)/2;
big_surf_idx = size(multilayer_struct.Vertices, 1)/2+1:size(multilayer_struct.Vertices, 1);

link_vectors = zeros(size(multilayer_struct.VertNormals));

for iVertex = 1:max(small_surf_idx)
    i_white = multilayer_struct.Vertices(small_surf_idx(iVertex), :);
    i_pial = multilayer_struct.Vertices(big_surf_idx(iVertex), :);
    
    link_vectors(iVertex, :) = (i_pial - i_white);%/norm(i_pial - i_white);
    link_vectors(iVertex + max(small_surf_idx), :) = (i_white - i_pial);%/norm(i_white - i_pial);
end

multilayer_struct.VertLinkV = link_vectors;