function res=QToomre(rr,vc,Md,Rd,alpha,zeta,GG)

%% find Toomre Q parameter for a galaxy embedded in a halo

[kappa,rnew]=epicyclic_frequency(rr,vc);

%units
%GG=G*km^-2*kpc^-1*Ms;

sigmaZ=alpha.*sqrt(zeta.*GG.*Md./Rd.*exp(-rnew./Rd));

sigmaS=Md./(2*pi*Rd.^2).*exp(-rnew./Rd);

res=sigmaZ.*kappa./(3.36.*GG.*sigmaS);



% 









%% what follows is the analytic expression for an isolated disk 
%J=bfunc(eta,1)+fg.*beta.^3.*bfunc(eta,beta)+2.*fb.*xi.^2./(eta.*(1+xi.*eta).^2);
%res=fac.*(pi/3.36).*eta.*sqrt(J).*epicyclic_frequency(eta,fg,beta,fb,xi)./(exp(-eta)+fg.*beta.^2.*exp(-beta.*eta));