function [result cosa] = Vr_bug(MPSec)

%%% Check for cached value
result = [];
%%%eval(sprintf('global Vr_RESULT_%d; result = Vr_RESULT_%d;',MPSec,MPSec));
%%%if (~isempty(result))
%%%    return;
%%%end

[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);
%Vxx = Vx(MPSec); Vyy = Vy(MPSec); Vzz = Vz(MPSec);

%%%% NOTE: Vyy is inversed!!!
[meshX, meshY, meshZ] = meshgrid(1:size(Vxx,1), size(Vyy,2):-1:1, 1:size(Vzz,3));
%convert to center origin coordinates
meshX = meshX - (size(Vxx,1)+1)/2;
meshY = meshY - (size(Vyy,2)+1)/2;
meshZ = meshZ - (size(Vzz,3)+1)/2;

norm = sqrt(meshX.^2 + meshY.^2 + meshZ.^2);
meshX = meshX ./ norm;
meshY = meshY ./ norm;
meshZ = meshZ ./ norm;

result = (Vxx.*meshX + Vyy.*meshY + Vzz.*meshZ);
V_abs = sqrt(Vxx.^2 + Vyy.^2 + Vzz.^2);
cosa = ((Vxx./V_abs).*meshX + (Vyy./V_abs).*meshY + (Vzz./V_abs).*meshZ);

%%% Cache the return value
%%%eval(sprintf('Vr_RESULT_%d = result;',MPSec));