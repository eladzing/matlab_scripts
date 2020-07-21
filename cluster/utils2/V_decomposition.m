function [Vr Vtheta Vphi] = V_decomposition(boxx,vcmflag,rvcm,hubflag,center)
% Calculates the (r,theta,phi) components of the velocity field based on
% V_Vcm()'s data.

if ~exist('center','var')
    center = [0,0,0];
end

if ~exist('rvcm','var')
    rvcm=[0 1];
end

if ~exist('hubflag','var')
    hubflag='hub';
end

%get velocities, after adding hubble flow & correcting for center of mass
%velocity
[Vxx Vyy Vzz] = get_velocities(boxx,vcmflag,rvcm,hubflag,center);

%if ~exist('cm')
%    cm = [0,0,0];
%end

%[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);

%calculate ^r
[meshY, meshX, meshZ] = meshgrid(1:size(Vyy,1), 1:size(Vxx,2), 1:size(Vzz,3));
%convert to center origin coordinates
meshX = meshX - (size(Vxx,1)+1)/2 -center(1);
meshY = meshY - (size(Vyy,2)+1)/2 -center(2);
meshZ = meshZ - (size(Vzz,3)+1)/2 -center(3);

meshR = sqrt(meshX.^2 + meshY.^2 + meshZ.^2);
% meshX = meshX ./ norm;
% meshY = meshY ./ norm;
% meshZ = meshZ ./ norm;

%calculate \theta and phi 
meshPhi = atan2(meshY,meshX);
meshTheta   = acos(meshZ./meshR);

Vr = (Vxx.*sin(meshTheta).*cos(meshPhi) + Vyy.*sin(meshTheta).*sin(meshPhi) + Vzz.*cos(meshTheta));
Vtheta = (Vxx.*cos(meshTheta).*cos(meshPhi) + Vyy.*cos(meshTheta).*sin(meshPhi) - Vzz.*sin(meshTheta));
Vphi = (-Vxx.*sin(meshPhi) + Vyy.*cos(meshPhi));
