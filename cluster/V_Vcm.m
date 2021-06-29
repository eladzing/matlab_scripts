function [Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec)

Vxx = Vx(MPSec);
Vyy = Vy(MPSec);
Vzz = Vz(MPSec);
Rhog = RHOG(MPSec);
M = sum(Rhog(:));

VcmX = sum(Vxx(:).*Rhog(:))/M;
VcmY = sum(Vyy(:).*Rhog(:))/M;
VcmZ = sum(Vzz(:).*Rhog(:))/M;

Vxx = Vxx - VcmX;
Vyy = Vyy - VcmY;
Vzz = Vzz - VcmZ;