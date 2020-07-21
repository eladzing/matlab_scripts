global cosmoStruct
global illUnits

% mcTNG=fofs.Group_M_Crit200.*illUnits.massUnit;
% mmTNG=fofs.Group_M_Mean200.*illUnits.massUnit;
% mask=mcTNG>1e12;
% 
% mcTNG=mcTNG(mask);
% mmTNG=mmTNG(mask);
% rcTNG=fofs.Group_R_Crit200(mask).*illUnits.lengthUnit;
% rmTNG=fofs.Group_R_Mean200(mask).*illUnits.lengthUnit;

% 
% lmcTNG=log10(mcTNG);
% lmmTNG=log10(mmTNG);


mh=10.^(12:0.05:15.3);
fgh=0.15;

rorefC=200*rho_crit(0,'cosmo',cosmoStruct)/1e9;
rorefM=200*rho_mean(0,'cosmo',cosmoStruct)/1e9;

rr=0.05:0.01:5;

for i=1:length(mh)
    
    cv=cvir_Mvir(mh(i),0);
    %crit
    
    
    hostC=NFW('mv',mh(i),'cc',cv,'fg',fgh,'crit');
    
    rom=meanDensity(hostC,rr);
    
    r2=interp1(rom,rr,rorefM);
    mm=mass(hostC,r2);
    
    crit2mean.Mc(i)=mh(i);
    crit2mean.Mm(i)=mm;
    crit2mean.Rc(i)=hostC.Rvir;
    crit2mean.Rm(i)=r2.*hostC.Rvir;
    crit2mean.cvc(i)=cv;
    crit2mean.cvm(i)=cvir_Mvir(mm,0);
    
    %mean
    
    
    hostM=NFW('mv',mh(i),'cc',cv,'fg',fgh,'mean');
    
    rom=meanDensity(hostM,rr);
    
    r2=interp1(rom,rr,rorefC);
    mm2=mass(hostM,r2);
    
    mean2crit.Mc(i)=mm2;
    mean2crit.Mm(i)=mh(i);
    mean2crit.Rm(i)=hostM.Rvir;
    mean2crit.Rc(i)=r2*hostM.Rvir;
    mean2crit.cvm(i)=cv;
    mean2crit.cvc(i)=cvir_Mvir(mm,0);
    
end
lmc1=log10(crit2mean.Mc);
lmm1=log10(crit2mean.Mm);

lmc2=log10(mean2crit.Mc);
lmm2=log10(mean2crit.Mm);




figure
%loglog(mcTNG(mask),mmTNG(mask),'.g')
%hold on
loglog(crit2mean.Mc,crit2mean.Mm,'r.')
hold on
loglog(mean2crit.Mc,mean2crit.Mm,'b.')

grid
loglog([1e12 1e15],[1e12 1e15],'--k')
xlabelmine('$M_\mathrm{crit}$');
ylabelmine('$M_\mathrm{mean}$');


figure
semilogx(mcTNG,mmTNG./mcTNG,'.g')
hold on
semilogx(crit2mean.Mc,crit2mean.Mm./crit2mean.Mc,'r.')

semilogx(mean2crit.Mc,mean2crit.Mm./mean2crit.Mc,'b.')

grid
%loglog([1e12 1e15],[1e12 1e15],'--k')
xlabelmine('$M_\mathrm{crit}$');
ylabelmine('$M_\mathrm{mean}/M_\mathrm{crit}$');

%% rvir

figure
loglog(rcTNG,rmTNG,'.g')
hold on
loglog(crit2mean.Rc,crit2mean.Rm,'r.')

loglog(mean2crit.Rc,mean2crit.Rm,'b.')

grid
loglog([180 2100],[180 2100],'--k')
xlim([180 2100])
ylim([180 2100])
xlabelmine('$R_\mathrm{crit}$');
ylabelmine('$R_\mathrm{mean}$');


figure
loglog(rcTNG,rmTNG,'.g')
hold on
loglog(crit2mean.Rc,crit2mean.Rm,'r.')

loglog(mean2crit.Rc,mean2crit.Rm,'b.')

grid
loglog([180 2100],[180 2100],'--k')
xlim([180 2100])
ylim([180 2100])
xlabelmine('$R_\mathrm{crit}$');
ylabelmine('$R_\mathrm{mean}$');




%% cv 

figure
%loglog(mcTNG(mask),mmTNG(mask),'.g')

plot(crit2mean.cvc,crit2mean.cvm,'r.')
hold on
plot(mean2crit.cvc,mean2crit.cvm,'b.')

grid
plot ([4 10],[4 10],'--k')
xlabelmine('$c_\mathrm{crit}$');
ylabelmine('$c_\mathrm{mean}$');




figure
loglog(crit2mean.Mc,crit2mean.Mm./crit2mean.Mc,'.')