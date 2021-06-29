function ff = calc_flux_sphere(MPSec,Vrcm,hubflag) 
%returns values in solarmass per year
%includes hubble flow
%Vrcm should include hubble flow.

Rhog = RHOG_sphere(MPSec);
Vrr  = Vr_sphere(MPSec);
ds   = ds_sphere(MPSec);

if ~exist('hubflag')
    hf=0;
else
    hf=(strcmpi(hubflag,'hub'));    
end

hubflo=hubble_sphere(MPSec);

if ~exist('Vrcm')
    Vrcm=zeros(size(Vrr));
end

Vrr = Vrr-Vrcm;
Vrr = Vrr+hf.*hubflo; 

% multiply by km / (1e6 * pc) * yr
ff    = Rhog.*Vrr.*ds.*(1.018120703e-12);

