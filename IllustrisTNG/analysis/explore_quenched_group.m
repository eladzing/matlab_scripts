
%x=[   10.1379     9.8269   9.9108   10.3255   11.6090   11.0018   10.0194   10.0194];
%y=[ 2.3425    1.3433    0.0921   -0.2410    2.1714    2.3425    2.3785    2.3785];

% x =[   11.1202   11.0018   10.5476   10.3057   10.0293   10.0737   10.1626 11.1202  ];
% y =[    2.1624    0.9922    0.3531    0.0291    0.3351    1.7034    2.3695 2.1624 ];

x=[11.1499   10.0046   10.0046   10.1872   11.1005   11.1252];
y=[2.2345    2.2975    0.2811   -0.1600    1.0822    2.2345];


% tc=tCoolStruct.(['in' gasField]).meanTcMW(1,tCoolStruct.galMask);
% gasTemp=tCoolStruct.(['in' gasField]).meanTempMW(1,tCoolStruct.galMask);
% 
% gasMach=tCoolStruct.(['in' gasField]).medianMach(1,tCoolStruct.galMask);
% gasMach2=tCoolStruct.(['in' gasField]).meanMachEW(1,tCoolStruct.galMask);
% gasEDisp=tCoolStruct.(['in' gasField]).EnergyDissipation(:,tCoolStruct.galMask);
% gasVDisp=tCoolStruct.(['in' gasField]).velDisp(:,tCoolStruct.galMask);

galMass=tCoolStruct.galMass;  % galaxy stellar mass
gasMass=tCoolStruct.inGal.gasMass(1,:);
gasEnt=tCoolStruct.inGal.meanEntMW(1,:);

%centralMask= subsInfo.isCentral(tCoolStruct.galMask);

sfr=subs.SubhaloSFRinRad;  % sfr in galaxy
ssfr=sfr./galMass; % + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));


%% get mask from polygon 
% define the mass range according to the selcted group (via polygon)
% find the rest of the popultion in the mass range , thus defining two
% distinct groups in the same mass range 

polyMask=inpolygon(log10(galMass),log10(gasEnt),x,y);
groupMask=polyMask & subsInfo.isCentral & tCoolStruct.galMask;
qGroupInd=find(groupMask);

massRange=[min(galMass(groupMask)) max(galMass(groupMask))];
massMask=galMass>= massRange(1) & galMass<=massRange(2);

oGroupMask=massMask & subsInfo.isCentral & tCoolStruct.galMask & ~groupMask ;
oGroupInd=find(oGroupMask);



%% show the distribution of host masses for the two groups 
qFofMass=fofs.Group_M_Crit200(subsInfo.hostFof(qGroupInd)+1).*illUnits.massUnit;
oFofMass=fofs.Group_M_Crit200(subsInfo.hostFof(oGroupInd)+1).*illUnits.massUnit;

nb=50;
binEdge=linspace(11,13.5,nb+1);
bins=binEdge(1:end-1)+0.5.*diff(binEdge);

qh=histcounts(log10(qFofMass),binEdge);
oh=histcounts(log10(oFofMass),binEdge);

figure
h(1)=stairs(bins,qh,'r','DisplayName','quenched','linewidth',2);
hold on
h(2)=stairs(bins,oh,'b','DisplayName','not-quenched','linewidth',2);
xlabelmine('$\log M_{\mathrm{200,c}}\,[\mathrm{M_\odot}]$')
ylabelmine('$N$')
grid
hl=legend(h);
set(hl,'interpreter','latex','fontsize',12,'location','NorthWest')

figure
h(1)=stairs(bins,qh./sum(qh),'r','DisplayName','quenched','linewidth',2);
hold on
h(2)=stairs(bins,oh./sum(oh),'b','DisplayName','not-quenched','linewidth',2);
xlabelmine('$\log M_{\mathrm{200,c}}\,[\mathrm{M_\odot}]$')
ylabelmine('$N / N_{\mathrm{total}}$')
grid
hl=legend(h);
set(hl,'interpreter','latex','fontsize',12,'location','NorthEast')

figure 
h(1)=stairs(bins,cumsum(qh)./sum(qh),'r','DisplayName','quenched','linewidth',2);
hold on
h(2)=stairs(bins,cumsum(oh)./sum(oh),'b','DisplayName','not-quenched','linewidth',2);
xlabelmine('$\log M_{\mathrm{200,c}}\,[\mathrm{M_\odot}]$')
ylabelmine('$N(<M)$')
grid
hl=legend(h);
set(hl,'interpreter','latex','fontsize',12,'location','NorthWest')


%% show the distribution of BH masses for the two groups 
qBHMass=subs.SubhaloMassInRadType(illustris.partTypeNum('bh')+1,qGroupInd).*illUnits.massUnit;
oBHMass=subs.SubhaloMassInRadType(illustris.partTypeNum('bh')+1,oGroupInd).*illUnits.massUnit;

nb=50;
binEdge=linspace(7,9.5,nb+1);
bins=binEdge(1:end-1)+0.5.*diff(binEdge);

qh=histcounts(log10(qBHMass),binEdge);
oh=histcounts(log10(oBHMass),binEdge);

figure
h(1)=stairs(bins,qh,'r','DisplayName','quenched','linewidth',2);
hold on
h(2)=stairs(bins,oh,'b','DisplayName','not-quenched','linewidth',2);
xlabelmine('$\log \Sigma M_{\mathrm{BH,Gal}}\,[\mathrm{M_\odot}]$')

ylabelmine('$N$')
grid
hl=legend(h);
set(hl,'interpreter','latex','fontsize',12,'location','NorthWest')

figure
h(1)=stairs(bins,qh./sum(qh),'r','DisplayName','quenched','linewidth',2);
hold on
h(2)=stairs(bins,oh./sum(oh),'b','DisplayName','not-quenched','linewidth',2);
xlabelmine('$\log \Sigma M_{\mathrm{BH,Gal}}\,[\mathrm{M_\odot}]$')
ylabelmine('$N / N_{\mathrm{total}}$')
grid
hl=legend(h);
set(hl,'interpreter','latex','fontsize',12,'location','NorthEast')

figure 
h(1)=stairs(bins,cumsum(qh)./sum(qh),'r','DisplayName','quenched','linewidth',2);
hold on
h(2)=stairs(bins,cumsum(oh)./sum(oh),'b','DisplayName','not-quenched','linewidth',2);
xlabelmine('$\log \Sigma M_{\mathrm{BH,Gal}}\,[\mathrm{M_\odot}]$')

ylabelmine('$N(<M)$')
grid
hl=legend(h);
set(hl,'interpreter','latex','fontsize',12,'location','NorthWest')


%% compare with mass-matched sample 
