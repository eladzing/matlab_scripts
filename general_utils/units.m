%% units in cgs 
Units.kb= 1.38066e-16; %erg/K
Units.muMass= 0.5926; %mass per particle for primordial composition
Units.mp= 1.672649e-24; % gram
Units.me= 9.10938356e-28; % gram 
Units.yr= 3.155815e7; %  sec 
Units.pc= 3.0856e18; % parsec
Units.km= 1e5; %center
Units.Ms= 1.9891e33; %gr
Units.xi= 0.5185185; % no. of electrons per particle - cosmology dependent
Units.G= 6.67e-8; % Gravitational Constant (cm^3 g^-1 s^-2)
Units.mm= Units.muMass*Units.mp;
Units.Mpc = Units.pc*1e6;
Units.kpc= Units.pc*1e3;
Units.ev= 1.602177e-12; % ergs
Units.Gyr=1e9*Units.yr;
Units.GG=Units.G*Units.km^-2*Units.kpc^-1*Units.Ms; % same as above in (kpc km/sec Msun^-1)
Units.cspeed= 2.99792458e10;  % speed of light in cm ;
Units.lightspeed=Units.cspeed;


Units.sigmaThomson=6.6524e-25; %Thomson cross sectionin cm^2  

%numerical factors converting to units of M_sun*(km/sec)^2 per Mpc^3
Units.factors.f_th=1.5.*Units.kb./(Units.muMass.*Units.mp)/Units.km.^2 ;%E_th factor converts to units cited above
Units.factors.f_pr=Units.kb./(Units.muMass.*Units.mp)/Units.km.^2 ;%E_th factor converts to units cited above
Units.factors.f_pv=Units.kb./(Units.muMass.*Units.mp.*Units.km.*Units.Mpc).*Units.yr; %units of M_sun*(km/sec)^2 per Mpc^3 per yr
Units.factors.fu=Units.km./(1e6.*Units.pc).*Units.yr;
Units.factors.f_cl=(Units.xi./(Units.muMass.*Units.mp)).^2*Units.yr*Units.Ms/(Units.Mpc)^3/Units.km^2*1.e-22;
Units.factors.f_eng=Units.G*Units.Ms/(Units.Mpc);

%other nuimerical factors: 
Units.factors.f_ent=(Units.kb/Units.ev/1e3)*(Units.Ms/Units.Mpc^3/(Units.muMass*Units.mp))^(-2/3); 
% converts entropy in simulation units (K, Msun/Mpc^3) (T/rho^(2/3) to keV*cm^2 

Units.header='All units in cgm. G in units of (cm^3 gr^-1 s^-2), GG units are (kpc (km / sec)^2 Msun^-1)';

