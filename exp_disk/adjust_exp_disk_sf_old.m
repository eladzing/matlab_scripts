function [sigma50,sigma5090]=adjust_exp_disk_sf_old(disc,rstrip)
% OBEOSLETE - see adjust_exp_disk_sf 
    
%%function to allow for an evolution of a galactic disc with gas
% an initial stellar component modeled as an exponential disk and a
% gasoeous component, also modeled as an exponential disk but truncated at
% a stripping radius rstrip (in kpc). 
% the function returns the new observed surface density measures. 

if ~isstruct(disc)
    error('adjust_exp_disk_sf: argument must be a structure')
end

rstrip=rstrip/disc.rd;
%% construct new mass profile
r=0:1e-4:10; % in units of the scale radius 
mstars=disc.mstar.*exp_disk_mass(r); % orginal stellar mass profile
gasmass=disc.mstar.*disc.fg.*exp_disk_mass(r,disc.beta); % original gas mass profile

if rstrip<r(end)
    ind=find(r>=rstrip);  
    gasmass(ind)=gasmass(ind(1));  % remove gas beyond rstrip 
end
totalgas=disc.mstar.*disc.fg.*exp_disk_mass(rstrip,disc.beta); % original gas mass profile

newmass=mstars+gasmass; % new mass profile 


% find new total mass 
bigRad=100;
if rstrip<bigRad
    fullmass=disc.mstar.*exp_disk_mass(bigRad)+totalgas;
else
    fullmass=disc.mstar.*exp_disk_mass(bigRad)+disc.mstar.*disc.fg.*exp_disk_mass(bigRad,disc.beta);
end

%% find r50 and r90
r50=interp1(newmass./fullmass,r,0.5).*disc.rd;
r90=interp1(newmass./fullmass,r,0.9).*disc.rd;

sigma50=0.5.*fullmass/(pi*r50^2);
sigma5090=0.4.*fullmass/(pi*(r90^2-r50^2));


