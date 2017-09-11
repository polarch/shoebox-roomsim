function echo2gridMap = get_echo2gridMap(echogram, grid_dirs_rad)
%QUANTISE_ECHOGRAM Quantises the reflections to specific rendering directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% QUANTISE_ECHOGRAM.M - 15/10/2011
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of reflections
nRefl = length(echogram.value);
% number of rendering directions
nGrid = size(grid_dirs_rad,1);

% unit vectors pointing to the rendering directions
[grid_xyz(:,1), grid_xyz(:,2), grid_xyz(:,3)] = sph2cart(grid_dirs_rad(:,1), grid_dirs_rad(:,2), 1);
% unit vectors pointing to the image sources
refl_coords = echogram.coords;
refl_coords = refl_coords./(sqrt(sum(refl_coords.^2,2))*ones(1,3));

for i=1:nRefl 
    [~, nearest] = min(sum((ones(nGrid,1)*refl_coords(i,:) - grid_xyz).^2, 2));
    echo2gridMap(i) = nearest;
end

end
