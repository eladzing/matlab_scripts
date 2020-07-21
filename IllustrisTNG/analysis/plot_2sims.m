%% construct plots with both TNG100 & TNG300
if ~skip
    cmap=brewermap(256,'Spectral');
    
    
    %% get tng 100
    fprintf(' *** get TNG100 data *** \n')
    
    
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/cooling_times_z0_TNG100.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_z0_TNG100.mat'])
    load([DEFAULT_MATFILE_DIR '/centralMask_TNG100.mat'])
    tcStr100=tCoolStruct;
    bhStr100=bhStruct;
    bp100=illustris.set_env('100');
    
    %% get data
    
    fprintf(' *** get TNG300 data *** \n')
    %global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/cooling_times_z0_TNG300.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_z0_TNG300.mat'])
    load([DEFAULT_MATFILE_DIR '/centralMask_TNG300.mat'])
    tcStr300=tCoolStruct;
    bhStr300=bhStruct;
    bp300=illustris.set_env('300');
    
    
    galMass100=tcStr100.galMass(centralMaskTNG100);
    galMass300=tcStr300.galMass(centralMaskTNG300);
    
    %gasMass100=tcStr100.(['in' gasField]).gasMass(:,centralMaskTNG100);
    %gasMass300=tcStr300.(['in' gasField]).gasMass(:,centralMaskTNG300);
    
    gasField='Gal';
    gasEnt100=tcStr100.(['in' gasField]).meanEntMW(1,centralMaskTNG100);
    gasEnt300=tcStr300.(['in' gasField]).meanEntMW(1,centralMaskTNG300);
    
    
    fprintf(' *** getting TNG100 SFR *** \n')
    subSFR100=illustris.groupcat.loadSubhalos(bp100,99,{'SubhaloSFRinRad'});
    sfr100=subSFR100(centralMaskTNG100);  % sfr in galaxy
    ssfr100=sfr100./galMass100 + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr100)));
    
    
    fprintf(' *** getting TNG300 SFR *** \n')
    subSFR300=illustris.groupcat.loadSubhalos(bp300,99,{'SubhaloSFRinRad'});
    sfr300=subSFR300(centralMaskTNG300);  % sfr in galaxy
    ssfr300=sfr300./galMass300 + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr300)));
    
    %% create contours 
    filt=fspecial('disk',8);
    xdata100=log10(galMass100);
    xdata300=log10(galMass300);
    xl=[9 13];
    
    ydata100=log10(gasEnt100);
    ydata300=log10(gasEnt300);
    yl=[-2.5 2.5];
    
    m100=ydata100>-Inf;
    m300=ydata300>-Inf;
    
    xdata100=xdata100(m100);
    ydata100=ydata100(m100);
    
    xdata300=xdata300(m300);
    ydata300=ydata300(m300);
    
    
    
    popCont100=plot_population_contour(xdata100,ydata100,'smooth',filt,'noplot');
    popCont300=plot_population_contour(xdata300,ydata300,'smooth',filt,'noplot');

    
    
    %% create color plot
    fprintf(' *** Creating data for plots *** \n')
    xdata=cat(2,xdata100,xdata300);
    %xl=[9 13];
       
    ydata=cat(2,ydata100,ydata300);
    %yl=[-2.5 2.5];
    cc=cat(2,log10(ssfr100(m100)),log10(ssfr300(m300)));
    cax=[min(cc) max(cc)];
    
    fprintf(' *** Building tree  *** \n')
    minLev=6;
    splitParam=100;
    tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
    celVal=points2tree(cc,tre,'median');
    
end

fprintf(' *** Plotting  *** \n')
ssfrLab='$\log \mathrm{sSFR},[\mathrm{yr^{-1}}] $';
% plot figure
figure

set(gcf,'position',[1432 421 1000 750],'Color','w')

ax1=axes;
plot2dTreeMap(celVal,tre,'cmap',(cmap),'minmax',cax,'fig',gcf);

hb=colorbar;barTitle(hb,ssfrLab,'fontsize',16)
set(hb,'fontsize',14);

colormap((cmap));caxis(cax);
xlim(xl);ylim(yl);
grid
ylabelmine('log Entropy')
xlabelmine('log Stellar Mass $\mathrm{M_\odot}$')

set(gca,'fontsize',14);

ax2=axes;
contour(popCont100.xx,popCont100.yy,popCont100.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
hold on
contour(popCont300.xx,popCont300.yy,popCont300.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];


printout_fig(gcf,'test100300','v')
