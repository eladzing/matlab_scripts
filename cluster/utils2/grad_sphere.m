function [gradR gradPhi gradTheta] = grad_sphere(MPsec, sphere_cube)
% calculates the gradient in spherical coordinates of the given spherical
% cube.

[mesh_r mesh_phi mesh_theta] = sphere_grid(MPsec);
mesh_r = single(mesh_r); mesh_phi = single(mesh_phi); mesh_theta = single(mesh_theta);

% d/dR
dR = mesh_r(2,1,1)-mesh_r(1,1,1);
gradR = (sphere_cube(3:end,:,:) - sphere_cube(1:end-2,:,:))/(2*dR);
gradR = gradR(:,2:end-1,2:end-1); %%% Chop to 254^3

% d/dPhi
dPhi = mesh_phi(:,2:end,:) - mesh_phi(:,1:end-1,:);
dPhi = dPhi(:,2:end,:) + dPhi(:,1:end-1,:);
gradPhi = ((sphere_cube(:,3:end,:) - sphere_cube(:,1:end-2,:))./dPhi)./mesh_r(:,2:end-1,:);
gradPhi = gradPhi(2:end-1,:,2:end-1); %%% Chop to 254^3

% d/dTheta
dTheta = mesh_theta(1,1,2) - mesh_theta(1,1,1);
gradTheta = ((sphere_cube(:,:,3:end) - sphere_cube(:,:,1:end-2))/(2*dTheta)) ./ (mesh_r(:,:,2:end-1) .* sin(mesh_phi(:,:,2:end-1)));
gradTheta = gradTheta(2:end-1,2:end-1,:); %%% Chop to 254^3

end