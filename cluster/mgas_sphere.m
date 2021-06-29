function mgas = mgas_sphere(MPSec) %returns values in solarmass

Rhog = RHOG_sphere(MPSec);
%%Vrr  = Vr_sphere(MPSec);
ds   = ds_sphere(MPSec);


mgas    = Rhog.*ds; 