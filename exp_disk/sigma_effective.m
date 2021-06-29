function sige=sigma_effective(Mdisk,rd)

%% calculate the effective sigma ( observational qunatity)
%Msat - mass in Msun
%rf - scale radius in kpc 

re=effective_rad(0.5)*rd; % half-mass radius

sige=0.5.*Mdisk./(pi.*re.^2);

