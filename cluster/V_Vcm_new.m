function [Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm_new(MPSec)

Vxx = Vx(MPSec);
Vyy = Vy(MPSec);
Vzz = Vz(MPSec);

% This loads VcmX VcmY and VcmZ
%load Vcm
global VCM
VcmX=VCM(1);
VcmY=VCM(2);
VcmZ=VCM(3);


Vxx = Vxx - VcmX;
Vyy = Vyy - VcmY;
Vzz = Vzz - VcmZ;