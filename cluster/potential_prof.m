function  [rpr pot]=potential_prof(rp)
% Calculates potential of cluster assuming spherical symmetry, based on
% total density profile for the current cluster (i.e. after new_env has 
% been invoked. function accepts as argument an array of radii.
% the function returns the the potential and an array of radii, after
% negative values have been removed.
% The potential is in units of (km/sec)^2

%cluster='CL101';
%type='csf';
%aexp='a1';


%% calculate potential profile for spherically symmetric density profile

%Gravitational constant with unit factors to set units to (km/sec)^2 
Gfac=-6.6726e-8*1.989e33/(1e6*3.0856e18)*1e-10; 

%read profiles from file
[~, rot] =read_RHO_Profiles(rp);
mt = read_MTOT_Profile(rp);

%get rid of negative values
rott=rot(rot>0 & mt>0);
rpr=rp(rot>0 & mt>0);
mtt=mt(rot>0 & mt>0);

%second integral 
rordr=4.*pi.*rpr.*rott; %integrand
% flip the array, preform cumulative integration and flip back
p1=fliplr(cumtrapz(rpr,fliplr(rordr))); 

pot=(mtt./rpr + p1).*Gfac; %first integral - M(<r)/r 


