function [ff, TT, SS] = flux(MPSec, R)

SPHERE_RES = 50*16; %*16;
INTERP = 'linear';
h = 0.8;

[Sx Sy Sz] = uni_sphere(SPHERE_RES);

%%%% Compensate Area
%%%n = SPHERE_RES;
%%%phi = (-n:2:n)'/n*pi/2;
%%%cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
%%%cosphi = repmat(cosphi, [1, n+1]);

%Sx = Sx(:); Sy = Sy(:); Sz = Sz(:);
RSx = R*Sx+(128.5);
RSy = R*Sy+(128.5);
RSz = R*Sz+(128.5);

%%% Optimized way to calculate Vrr. Probably works, but not sure.
%Vxx  = interp3(Vx(MPSec),   RSx, RSy, RSz, INTERP);
%Vyy  = interp3(Vy(MPSec),   RSx, RSy, RSz, INTERP); 
%Vzz  = interp3(Vz(MPSec),   RSx, RSy, RSz, INTERP); 
%Vrr = dot([Vxx Vyy Vzz],[Sx (-1)*Sy Sz],2);

Rhog = interp3(RHOG(MPSec), RSx, RSy, RSz, INTERP);
Vrr  = interp3(Vr(MPSec),   RSx, RSy, RSz, INTERP);
ff    = Rhog.*Vrr*((R*(MPSec/h)/256)^2); %%%.*cosphi;

TT    = interp3(T(MPSec),   RSx, RSy, RSz, INTERP);
SS    = interp3(S(MPSec),   RSx, RSy, RSz, INTERP);

%ff = ff(:);
%TT = TT(:);
%SS = SS(:);