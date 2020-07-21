function [sigma50,sigma5090]=adjust_exp_disk_sf(rstrip,fg,beta)
%% function to allow for an evolution of a galactic disc with gas
% an initial stellar component modeled as an exponential disk and a
% gasoeous component, also modeled as an exponential disk but truncated at
% a stripping radius rstrip (in units of stellar scale radius). 
% the function returns the factors to multiply Sigma_s (Ms/(2*pi*rd^2)) to obtain the 
% new observed surface density measures

%% construct new mass profile
r=0:1e-4:10; % in units of the scale radius 
mstars=exp_disk_mass(r); % orginal stellar mass profile
gasmass=fg.*exp_disk_mass(r,beta); % original gas mass profile

if rstrip<r(end)
    ind=find(r>=rstrip);  
    gasmass(ind)=gasmass(ind(1));  % remove gas beyond rstrip 
end
totalgas=fg.*exp_disk_mass(rstrip,beta); % original gas mass profile

newmass=mstars+gasmass; % new mass profile 

% find new total mass 
bigRad=100;
if rstrip<bigRad
    fullmass=1+totalgas;
else
    fullmass=1+fg.*exp_disk_mass(bigRad,beta);
end

%% find r50 and r90
r50=interp1(newmass./fullmass,r,0.5);
r90=interp1(newmass./fullmass,r,0.9);

sigma50=0.5.*fullmass/(pi*r50^2)*2*pi;
sigma5090=0.4.*fullmass/(pi*(r90^2-r50^2))*2*pi;


