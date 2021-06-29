function warped = warp_mat3d(cube, transform_matrix)
%warp_rot - warps an imgage using back-warping. The homogeneous transformation 
% matrix is given, as a 3x3 matrix.
%
%cube - grayscale image to warp
%transform_matrix - the transformation matrix

%prepare a mesh grid of lookups into the original image
[meshX, meshY, meshZ] = meshgrid(1:size(cube,1), 1:size(cube,2), (size(cube,3)+1)/2);
%convert to center origin coordinates
meshX = meshX - (size(cube,1)+1)/2;
meshY = meshY - (size(cube,2)+1)/2;
meshZ = meshZ - (size(cube,3)+1)/2;

%create a homogeneous coordinates representation of the grid.
homoTarget = [transpose(meshX(:)); transpose(meshY(:)); transpose(meshZ(:)); ...
    repmat([1],1,size(meshY(:)))];

%apply the transformation matrix on the coordinates.
solutions = transform_matrix * homoTarget;

%convert solutions back to left-corner origin non-homogeneous coordinates.
solutions(1,:) = solutions(1,:) + (size(cube,1)+1)/2;
solutions(2,:) = solutions(2,:) + (size(cube,2)+1)/2;
solutions(3,:) = solutions(3,:) + (size(cube,3)+1)/2;

%convert homogeneous coordinates to X, Y matrices
XI = reshape(solutions(1,:), [size(cube,1) size(cube,2)]);
YI = reshape(solutions(2,:), [size(cube,1) size(cube,2)]);
ZI = reshape(solutions(3,:), [size(cube,1) size(cube,2)]);

%interpolate (bilinear interpolation)
warped = interp3(cube, XI, YI, ZI);

%set all invalid values (pixels that we don't have image data for them) to 
% 0 (black).
%warped(isnan(warped)) = 0;