function sfr=sfr_proxy_sigma(mgas,rval)%,beta) %,rmax)

%rmax=5.*rd;

A=2.5e-4; % in units of solar mass per yer per kpc^-2

% sigma gas is in solar mass to pc^2
%siggas=mgas.*exp_disk_mass(rmax./rd,beta)./(pi.*(1e3.*rmax).^2); 
siggas=mgas./(pi.*(1e3.*rval).^2); 

sfr=A.*(siggas).^1.4;

%ssfr=sfr./ms.*exp_disk_mass(rmax./rd,1));
%ssfr=sfr./ms;

