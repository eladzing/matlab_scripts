function [result] = Vr_full(MPSec)
% Calculates the radial component of the velocity based on V_Vcm()'s data.

global VCM;
%global zred;
%global hub;

cm=[0 0 0];
% get hubble flow 
[hubX hubY hubZ] = hubble_flow(MPSec); 
hf=1.0;
Vxx = Vx(MPSec)+hf.*hubX-VCM(1);
Vyy = Vy(MPSec)+hf.*hubY-VCM(2);
Vzz = Vz(MPSec)+hf.*hubZ-VCM(3);

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
