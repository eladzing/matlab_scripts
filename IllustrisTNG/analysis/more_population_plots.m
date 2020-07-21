ssfr=illustris.utils.calc_ssfr(subs);


ydata=log10(gasEnt(galaxyMask));
yl=[-2 3];
cc=log10(ssfr(galaxyMask));
cax=[min(cc) max(cc)];


%% ratio 
bhrat=log10(bhRM(galaxyMask)./bhQM(galaxyMask));
xdata=bhrat;xl=[-7 0];
minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=log10(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
ylabelmine('log Entropy')
xlabelmine('log (KM/QM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);


%% QM

bhrat=log10(bhQM(galaxyMask));
xdata=bhrat;
xl=[5 9];
minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=log10(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',flipud(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(flipud(cmap));caxis(cax);
grid
ylabelmine('log Entropy')
xlabelmine('log (QM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);


%% KM
mask=galaxyMask & bhRM>0;
ydata=log10(gasEnt(mask));
yl=[-2.5 2.5];
bhrat=log10(bhRM(mask));
xdata=bhrat;
xl=[0 9];

cc=log10(bhQM(mask));
cax=[5 9];

minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
%celVal=log10(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(cmap);caxis(cax);
grid
ylabelmine('log Entropy')
xlabelmine('log (KM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,'QM','fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% bh mass 

%cc=(ssfr(galaxyMask));
ydata=log10(gasEnt(gm));
yl=[-2.5 2.5];
cc=(ssfr(gm));

bhrat=log10(bhm(gm).*illUnits.massUnit);
xdata=bhrat;
xl=[5 10.1];
minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=log10(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',flipud(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(flipud(cmap));caxis(cax);
grid
ylabelmine('log Entropy')
xlabelmine('log $M_\mathrm{BH}$')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% QM vs. RM 
ydata=log10(bhRM(galaxyMask));
yl=[-1 9];


xdata=log10(bhQM(galaxyMask));
xl=[6 9];

cc=log10(ssfr(galaxyMask));
cax=[min(cc) max(cc)];

minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',flipud(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(flipud(cmap));caxis(cax);
grid
ylabelmine('log RM')
xlabelmine('log (QM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% M vs. RM 
xdata=log10(galMass(galaxyMask));

xl=[9 13];


ydata=log10(bhRM(galaxyMask));
yl=[5 9];

cc=log10(ssfr(galaxyMask));
cax=[min(cc) max(cc)];

minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',flipud(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(flipud(cmap));caxis(cax);
grid
xlabelmine('log M')
ylabelmine('log (RM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% M vs. QM 
xdata=log10(galMass(galaxyMask));

xl=[9 13];


ydata=log10(bhQM(galaxyMask));
yl=[5 9];

cc=log10(ssfr(galaxyMask));
cax=[min(cc) max(cc)];

minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',flipud(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(flipud(cmap));caxis(cax);
grid
xlabelmine('log M')
ylabelmine('log (QM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% Mb vs. RM 
xdata=log10(bhMass(galaxyMask));

xl=[7 11];


ydata=log10(bhRM(galaxyMask));
yl=[5 9];

cc=log10(ssfr(galaxyMask));
cax=[min(cc) max(cc)];

minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',flipud(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(flipud(cmap));caxis(cax);
grid
xlabelmine('log Mb')
ylabelmine('log (RM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% M vs. QM 
xdata=log10(bhMass(galaxyMask));

xl=[7 11];


ydata=log10(bhQM(galaxyMask));
yl=[5 9];

cc=log10(ssfr(galaxyMask));
cax=[min(cc) max(cc)];

minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',flipud(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(flipud(cmap));caxis(cax);
grid
xlabelmine('log Mb')
ylabelmine('log (QM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);

%% Mb vs. RM 
xdata=log10(bhMass(galaxyMask));

xl=[7 11];


ydata=log10(bhRM(galaxyMask)./bhQM(galaxyMask));
yl=[-6 0];

cc=log10(ssfr(galaxyMask));
cax=[min(cc) max(cc)];

minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',flipud(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap(flipud(cmap));caxis(cax);
grid
xlabelmine('log Mb')
ylabelmine('log (RM/QM)')
ssfrLab='$\log \mathrm{sSFR}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);



%%   density plots 

xdata=log10(galMass(galaxyMask));
xl=[9 13];

ydata=log10(gasEnt(galaxyMask));
yl=[-2 3];

%cc=log10(ssfr(galaxyMask));
cc=log10(medianRho(galaxyMask));
cax=[min(cc) max(cc)];

minLev=6;
splitParam=20;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
celVal=points2tree(cc,tre,'median');
celVal=(celVal);

figure
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);
hb=colorbar;

colormap((cmap));caxis(cax);
grid
xlabelmine('log M')
ylabelmine('log K')
ssfrLab='$\log \mathrm{n}_{\mathrm{Gal}},[\mathrm{yr^{-1}}] $';
barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);
set(gca,'fontsize',14);



