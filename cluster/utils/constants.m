%% general constants used often, in c.g.s unless noted otherwise


kb=1.38066e-16; %k_boltzman :: erg/K
mu=0.5926; %mass per particle for primordial composition
mp=1.672649e-24; % mass of protons in grams
yr=3.155815e7; %  1 year in sec 
pc=3.0856e18; % 1 paraec in cm
km=1e5; %1 kilometer in cm
Ms=1.989e33; %1 solar mass in gr
xi=0.5185185; % ionization fraction (?) used in thermal calculations 

%% numerical factors converting energy types 
% from erg/cm^3  to units of M_sun*(km/sec)^2 per Mpc^3

f_th=1.5.*kb./(mu.*mp)/km.^2 ;%E_th factor converts to units cited above
f_pv=kb./(mu.*mp.*km.*1.0e6.*pc).*yr; %units of M_sun*(km/sec)^2 per Mpc^3 per yr
fu=km./(1e6.*pc).*yr; % converts velocity 
f_cl=(xi./(mu.*mp)).^2*yr*Ms/(1e6*pc)^3/km^2*1.e-22; % cooling constants and conversion 
