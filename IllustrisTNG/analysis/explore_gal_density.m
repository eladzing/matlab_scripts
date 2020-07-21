%% explore density plots 

% general things 
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
mbhLab='$\log M_\mathrm{BH}\,[\mathrm{M_\odot}]$';
mhostLab='$\log M_\mathrm{200}\,[\mathrm{M_\odot}]$';
minLev=6;
splitParam=20;
cmap=brewermap(256,'*Spectral');


%% data stuff 
gasEnt=tCoolStruct.inGal.meanEntMW(1,:);
meanRho=tCoolStruct.inGal.meanDensN(1,:);
avgRho=tCoolStruct.inGal.avgDensN(1,:);
galMass=tCoolStruct.galMass;
bhMass=subs.SubhaloMassInRadType(illustris.partTypeNum('bh')+1,:).*illUnits.massUnit;
bhQM=bhStruct.inGal.cumEngQM;
bhRM=bhStruct.inGal.cumEngRM;
hostMass=fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit;

ssfr=illustris.utils.calc_ssfr(subs);
% masks 
mask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',10^9.5,'centrals');

%% Galaxy

galaxyMask=mask & tCoolStruct.inGal.gasMass(1,:)>0;

entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$','Gal');

%% mean rho 

xdata=log10(galMass(galaxyMask));
xl=[9 13];

ydata=log10(gasEnt(galaxyMask));
yl=[-3 3];

%cc=log10(ssfr(galaxyMask));
cc=log10(meanRho(galaxyMask));
%cax=[-4 0];
cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
xlabelmine(mstarLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{Gal,mean}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% avg rho 


%cc=log10(ssfr(galaxyMask));
cc=log10(avgRho(galaxyMask));
cax=[-4 0];
%cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
xlabelmine(mstarLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{Gal,avg}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);


%% vs Mbh 

%% mean rho 

xdata=log10(bhMass(galaxyMask));
xl=[6 10.1];

ydata=log10(gasEnt(galaxyMask));
yl=[-3 3];

%cc=log10(ssfr(galaxyMask));
cc=log10(meanRho(galaxyMask));
cax=[-4 0];
%cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
xlabelmine(mbhLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{Gal,mean}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% avg rho 


%cc=log10(ssfr(galaxyMask));
cc=log10(avgRho(galaxyMask));
cax=[-4 0];
%cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
xlabelmine(mbhLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{Gal,avg}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);







%% ====================================
 
%% CGM 

gasEnt=tCoolStruct.inCGM.meanEntMW(1,:);
meanRho=tCoolStruct.inCGM.meanDensN(1,:);
avgRho=tCoolStruct.inCGM.avgDensN(1,:);

galaxyMask=mask & tCoolStruct.inCGM.gasMass(1,:)>0;
entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$','CGM');


%% mean rho 

xdata=log10(galMass(galaxyMask));
xl=[9 12.5];

ydata=log10(gasEnt(galaxyMask));
yl=[-3 2.5];

%cc=log10(ssfr(galaxyMask));
cc=log10(meanRho(galaxyMask));
cax=[-5 -1];
%cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;
ylim(yl)
colormap((cmap));caxis(cax);
grid
xlabelmine(mstarLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{CGM,mean}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% avg rho 


%cc=log10(ssfr(galaxyMask));
cc=log10(avgRho(galaxyMask));
cax=[-5 -1];
%cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;
ylim(yl)
colormap((cmap));caxis(cax);
grid
xlabelmine(mstarLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{CGM,avg}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);


%% vs Mbh 

%% mean rho 

xdata=log10(bhMass(galaxyMask));
xl=[6 10.1];

ydata=log10(gasEnt(galaxyMask));
yl=[-3 3];

%cc=log10(ssfr(galaxyMask));
cc=log10(meanRho(galaxyMask));
cax=[-5 -1];
%cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
xlabelmine(mbhLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{CGM,mean}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% avg rho 

%cc=log10(ssfr(galaxyMask));
cc=log10(avgRho(galaxyMask));
cax=[-5 -1];
%cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
xlabelmine(mbhLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{CGM,avg}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% by halo mass 

%% avg rho 
xdata=log10(hostMass(galaxyMask));
xl=[9 15];

%cc=log10(ssfr(galaxyMask));
cc=log10(avgRho(galaxyMask));
cax=[-5 -1];
%cax=[min(cc) max(cc)];


tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
xlabelmine(mhostLab) 
ylabelmine(entLab)
ssfrLab='$\log \mathrm{n}_{\mathrm{CGM,avg}},[\mathrm{cm^{-3}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);




%% vs  bh ratio 



bhrat=log10(bhRM(galaxyMask)./bhQM(galaxyMask));
xdata=bhrat;xl=[-7 0];
