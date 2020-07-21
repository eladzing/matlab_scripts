function [result] = Vr(MPSec, cm)
% Calculates the radial component of the velocity based on V_Vcm()'s data.

if ~exist('cm')
    cm = [0,0,0];
end

Vxx = Vx(MPSec);
Vyy = Vy(MPSec);
Vzz = Vz(MPSec);

[meshY, meshX, meshZ] = meshgrid(1:size(Vyy,1), 1:size(Vxx,2), 1:size(Vzz,3));
%convert to center origin coordinates
meshX = meshX - (size(Vxx,1)+1)/2 -cm(1);
meshY = meshY - (size(Vyy,2)+1)/2 -cm(2);
meshZ = meshZ - (size(Vzz,3)+1)/2 -cm(3);

norm = sqrt(meshX.^2 + meshY.^2 + meshZ.^2);
meshX = meshX ./ norm;
meshY = meshY ./ norm;
meshZ = meshZ ./ norm;

result = (Vxx.*meshX + Vyy.*meshY + Vzz.*meshZ);
