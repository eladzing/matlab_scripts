function ff = energy_sphere(MPSec) %returns values in solarmass * (km/sec)^2 per year

Rhog = RHOG_sphere(MPSec);
Vrr  = Vr_sphere(MPSec);
ds   = ds_sphere(MPSec);

% multiply by km / (1e6 * Mpc) * yr
ff    = 0.5 .* Rhog.*(Vrr)^2.*ds.*(1.018120703e-12);