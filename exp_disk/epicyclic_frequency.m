function [kappa,rNew]=epicyclic_frequency(rr,vc)
%% epicyclic frequency given the circular velocity (in units of km/sec) and radial profile (kpc)
% since we are using employing a numerical derivative we must reset the
% profile to excluse the first and last point. 
units

rKm=rr; 

[dvdr,rKm2]=derive1(vc,rKm);

indx=2:length(rr)-1;

v=vc(indx);
rNew=rr(indx);

kappa=sqrt(2.*(v.^2./rKm2.^2+(v./rKm2).*dvdr)) ; %in units of (km/sec) kpc^-1





%% what follows is the analytical expression for the case of an isolated exponential disk 

%J=bfunc(eta,1)+fg.*beta.^3.*bfunc(eta,beta)+2.*fb.*xi.^2./(eta.*(1+xi.*eta).^2);

%dJ=bfunc_derivative(eta,1)...
%    +fg.*beta.^3.*bfunc_derivative(eta,beta)...
%    -2.*fb.*xi.^2.*(1+3.*xi.*eta)./(eta.^2.*(1+xi.*eta).^3);

%res=sqrt(4.*J+eta.*dJ);