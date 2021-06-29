function [Vr Vtheta Vphi] = V_decomposition(MPSec, cm)
% Calculates the (r,theta,phi) components of the velocity field based on
% V_Vcm()'s data.

if ~exist('cm')
    cm = [0,0,0];
end

[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);

%calculate ^r
[meshY, meshX, meshZ] = meshgrid(1:size(Vyy,1), 1:size(Vxx,2), 1:size(Vzz,3));
%convert to center origin coordinates
meshX = meshX - (size(Vxx,1)+1)/2 -cm(1);
meshY = meshY - (size(Vyy,2)+1)/2 -cm(2);
meshZ = meshZ - (size(Vzz,3)+1)/2 -cm(3);

meshR = sqrt(meshX.^2 + meshY.^2 + meshZ.^2);
% meshX = meshX ./ norm;
% meshY = meshY ./ norm;
% meshZ = meshZ ./ norm;

%calculate \theta
meshTheta = atan2(meshY,meshX);
meshPhi   = acos(meshZ./meshR);

Vr = (Vxx.*cos(meshTheta).*sin(meshPhi) + Vyy.*sin(meshTheta).*sin(meshPhi) + Vzz.*cos(meshPhi));
Vtheta = (-Vxx.*sin(meshTheta) + Vyy.*cos(meshTheta));
Vphi = (Vxx.*cos(meshTheta).*cos(meshPhi) + Vyy.*sin(meshTheta).*cos(meshPhi) - Vzz.*sin(meshPhi));
