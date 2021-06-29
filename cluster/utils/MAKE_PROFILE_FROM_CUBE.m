function result = MAKE_PROFILE_FROM_CUBE(data)
% Returns a sum profile for cartesian cubes of extensive parameters.  
% The profile is the sum of all values at a given radius, and not the average. 
% The calculation is done in a manner that preserves the sum of all data 
% within a sphere, Unlike using using cart2sphere and then summing, which 
% resamples and does not preserve cumulative values. Useful when the total 
% 'mass' should be preserved, or when calculating integrals.
%
% @param data    The cartesian data cube.
%
% @returns     The sum profile
%
narginchk(1,1);

% Create meshR
[meshY, meshX, meshZ] = meshgrid(1:size(data,1), 1:size(data,2), 1:size(data,3));
%convert to center origin coordinates
meshX = meshX - (size(data,1)+1)/2;
meshY = meshY - (size(data,2)+1)/2;
meshZ = meshZ - (size(data,3)+1)/2;
meshR = sqrt(meshX.^2 + meshY.^2 + meshZ.^2);
clear meshX meshY meshZ

lowerCell     = int32(floor(meshR * 2));
upperFraction = (meshR * 2) - floor(meshR * 2);

result = zeros([1 257]);
data = double(data);
for i = 1:256
    mask = find(lowerCell == i);
    result(i)   = result(i)   + sum(data(mask) .* (1 - upperFraction(mask)));
    result(i+1) = result(i+1) + sum(data(mask) .* upperFraction(mask));
end

result = result(1:256);