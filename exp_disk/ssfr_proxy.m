function ssfr=ssfr_proxy(ms,mgas,rmax)%beta,rd) %,rmax)

%rmax=5.*rd;

A=(2.5).*1e-4;

% sigma gas is in solar mass to pc^2
%siggas=mgas.*exp_disk_mass(rmax./rd,beta)./(pi.*(1e3.*rmax).^2); 
siggas=mgas./(pi.*(1e3.*rmax).^2); 

sfr=(pi.*(rmax).^2).*A.*(siggas).^(1.4+0.15);


%ssfr=sfr./ms.*exp_disk_mass(rmax./rd,1));
ssfr=sfr./ms;

