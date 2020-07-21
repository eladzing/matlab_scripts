function sfr=sfr_proxy(mgas,rd,beta,rmax)


A=2.5e-4;

% sigma gas is in solar mass to pc^2
siggas=mgas.*exp_disk_mass(rmax./rd,beta)./(pi.*(1e3.*rmax).^2); 

sfr=pi.*rmax.^2.*A.*(siggas).^1.4;


