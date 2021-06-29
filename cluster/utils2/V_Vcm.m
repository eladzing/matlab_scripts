function [Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec,hubflag,cm)
% Calculates the center of mass velocity and returns the corrected velocity
% field. if hubflag is present, add hubble flow
if ~exist('hubflag')
    hf=0;
else
    hf=(strcmpi(hubflag,'hub'));    
end

if ~exist('cm','var')
    cm = [0,0,0];
end

[hubX hubY hubZ] = hubble_flow(MPSec,cm); 

Vxx = Vx(MPSec)+hf.*hubX;
Vyy = Vy(MPSec)+hf.*hubY;
Vzz = Vz(MPSec)+hf.*hubZ;

Rhog = RHOG(MPSec);
M = sum(Rhog(:));

VcmX = sum(Vxx(:).*Rhog(:))/M;
VcmY = sum(Vyy(:).*Rhog(:))/M;
VcmZ = sum(Vzz(:).*Rhog(:))/M;

Vxx = Vxx - VcmX;
Vyy = Vyy - VcmY;
Vzz = Vzz - VcmZ;