function [rvir mvir vvir tvir deltav]=calc_virials(delv)
%% Calculate Virial qunatities for different values of delta_vir.  


%% units in cgs 
kb= 1.38066e-16; %erg/K
mu= 0.5926; %mass per particle for primordial composition
mp= 1.672649e-24; % gram
pc= 3.0856e18; % parsec
km= 1e5; %center
Ms= 1.9891e33; %gr
G= 6.67e-8; % Gravitational Constant (cm^3 g^-1 s-2)

global zred

if ~exist('delv','var')
    deltav=deltavir(zred);
else
    deltav=delv;
end


[rvir mvir]=get_virs_from_profile(deltav);

vvir=sqrt(G*(mvir*Ms)/(rvir.*1e6*pc))/km; %in km/sec
tvir=(vvir*km)^2*mu*mp/(2*kb); % inKelvin



    