function warped = warp_mat3d_fast(cube, transform_matrix, FACTOR, DEPTH)
%warp_rot - warps an imgage using back-warping. The homogeneous transformation 
% matrix is given, as a 3x3 matrix.
%
%cube - grayscale image to warp
%transform_matrix - the transformation matrix
STEP = (size(cube,1) - 1) / (size(cube,1)*FACTOR - 1); %1/FACTOR;

%prepare a mesh grid of lookups into the original image
[meshX, meshY, meshZ] = meshgrid(1:STEP:size(cube,1), 1:STEP:size(cube,2), DEPTH);%(size(cube,3)+1)/2);
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
%XI = reshape(solutions(1,:), [size(cube,1)*FACTOR size(cube,2)*FACTOR]);
%YI = reshape(solutions(2,:), [size(cube,1)*FACTOR size(cube,2)*FACTOR]);
%ZI = reshape(solutions(3,:), [size(cube,1)*FACTOR size(cube,2)*FACTOR]);
%warped = interp3(cube, XI, YI, ZI);

%interpolate (bilinear interpolation)
%[x y z] = meshgrid(1:size(cube,1),1:size(cube,2),1:size(cube,3));
global xmesh;
global ymesh;
global zmesh;
warped = interp3(xmesh,ymesh,zmesh,cube, solutions(1,:), solutions(2,:), solutions(3,:),'nearest'); 
warped = reshape(warped, [size(cube,1)*FACTOR size(cube,2)*FACTOR]);

%set all invalid values (pixels that we don't have image data for them) to 
% 0 (black).
warped(isnan(warped)) = 0;