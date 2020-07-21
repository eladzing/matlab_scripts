function qq = flux_sphere(MPSec) %returns values in solarmass per year

Rhog = RHOG_sphere(MPSec);
T6  = T_sphere(MPSec)./1e6;
%ds   = ds_sphere(MPSec);

zmet=0.3;

lmb23=6.0.*(zmet/0.3)^0.7.*T6.^(-1)+0.2.*T6.^0.5;
qq=1.8991524.425.*Rhog.*lmb23; %in units of erg per sec per gr
