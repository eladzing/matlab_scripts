%% fix Jeans tot
Jean=JeanTot_a1;
units;
%rp=0.01:0.001:4/0.7; %till size of biggest box
%rp=rp(2:end-1);
for j=1:10
    
    %rp=Jean(j).rp;
    
    sig=Jean(j).sigma;
    rog=Jean(j).rog;
    
    sigR0=cat(2,sig(1).sigSqR(1:end-1),sig(2).sigSqR(128:end-1),sig(3).sigSqR(128:end-1),sig(4).sigSqR(128:end));
    sigT0=cat(2,sig(1).sigSqT(1:end-1),sig(2).sigSqT(128:end-1),sig(3).sigSqT(128:end-1),sig(4).sigSqT(128:end));
    sigP0=cat(2,sig(1).sigSqP(1:end-1),sig(2).sigSqP(128:end-1),sig(3).sigSqP(128:end-1),sig(4).sigSqP(128:end));
    
    vrP0=cat(2,sig(1).vrs(1:end-1),sig(2).vrs(128:end-1),sig(3).vrs(128:end-1),sig(4).vrs(128:end));
    vthP0=cat(2,sig(1).vths(1:end-1),sig(2).vths(128:end-1),sig(3).vths(128:end-1),sig(4).vths(128:end));
    vphP0=cat(2,sig(1).vphs(1:end-1),sig(2).vphs(128:end-1),sig(3).vphs(128:end-1),sig(4).vphs(128:end));
    
    rpSig0=cat(2,sig(1).rp(1:end-1),sig(2).rp(128:end-1),sig(3).rp(128:end-1),sig(4).rp(128:end));
    
    
    sigR=interp1(rpSig0,sigR0,rp);
    sigT=interp1(rpSig0,sigT0,rp);
    sigP=interp1(rpSig0,sigP0,rp);
    
    vrP=interp1(rpSig0,vrP0,rp);
    vthP=interp1(rpSig0,vthP0,rp);
    vphP=interp1(rpSig0,vphP0,rp);
    
    fac=(rog.*Ms./Mpc^3)./(rp.*Mpc);   
    isotrop=fac.*(2.*sigR-sigT-sigP);
    spheric=fac.*(vthP.^2+vphP.^2);
    
    
    
    
    Jean(j).sigR=sigR;
    Jean(j).sigT=sigT;
    Jean(j).sigP=sigP;
    Jean(j).vrP=vrP;
    Jean(j).vthP=vthP;
    Jean(j).vphP=vphP;
    Jean(j).isotrop=isotrop;
    Jean(j).spheric=spheric;
end
