%% Plot results for gas properties in Centrals in TNG
%% load data
global simDisplayName

if readFlag
    
    global DEFAULT_MATFILE_DIR
    %load([DEFAULT_MATFILE_DIR '/cooling_times_z0_' simDisplayName '.mat'])
    load([DEFAULT_MATFILE_DIR '/gasProperties_snp' num2str(snap) '_' simDisplayName '.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_snp' num2str(snap) '_' simDisplayName '.mat'])
    %     bhStruct.inGal=illustris.utils.read_catalog('bh_inGal','folder','bhProps');
    %     tCoolStruct.inGal=illustris.utils.read_catalog('gasProps_inGal','folder','gasProperties');
    
    loadFofSub
    %
    %         global DRACOFLAG
    %         if DRACOFLAG
    %             fofs=illustris.groupcat.loadHalos(bp,snap);
    %             subs=illustris.groupcat.loadSubhalos(bp,snap);
    %
    %         else
    %             loadFofSubTNG100
    %         end
    %
    %     fofs=illustris.groupcat.loadHalos(bp,snap);
    %     subs=illustris.groupcat.loadSubhalos(bp,snap);
    
end
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
fofs=illustris.utils.addTvirFofs( fofs);
units; % load general unit structure in cgs.
illustris.utils.set_illUnits(snap)
global illUnits

%% define what we are plotting
% identify centrals
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

%gasField='Gal';

%global simDisplayName

% get usefel stuff
galMass=tCoolStruct.galMass(tCoolStruct.galMask);  % galaxy stellar mass
gasMass=tCoolStruct.(['in' gasField]).gasMass(:,tCoolStruct.galMask);
sfr=subs.SubhaloSFRinRad(tCoolStruct.galMask);  % sfr in galaxy
ssfr=illustris.utils.calc_ssfr(subs);% sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));
ssfr=ssfr(tCoolStruct.galMask);

switch snap
    case 99
        ssfrThresh=1e-11;
    case 67
        ssfrThresh=10^(-10.3);
    otherwise
        error('sSFR threshold for quenching not defined for this snapshot');
end
qMask=ssfr<=ssfrThresh;

% sfre(1,:)=sfr./gasMass(1,:);
% sfre(2,:)=sfr./gasMass(1,:);



tc=tCoolStruct.(['in' gasField]).meanTcMW(1,tCoolStruct.galMask);
gasTemp=tCoolStruct.(['in' gasField]).meanTempMW(1,tCoolStruct.galMask);
gasEnt=tCoolStruct.(['in' gasField]).meanEntMW(1,tCoolStruct.galMask);
% gasMach=tCoolStruct.(['in' gasField]).medianMach(1,tCoolStruct.galMask);
% gasMach2=tCoolStruct.(['in' gasField]).meanMachEW(1,tCoolStruct.galMask);
% gasEDisp=tCoolStruct.(['in' gasField]).EnergyDissipation(:,tCoolStruct.galMask);
% gasVDisp=tCoolStruct.(['in' gasField]).velDisp(:,tCoolStruct.galMask);
gasDens=tCoolStruct.(['in' gasField]).meanDensN(1,tCoolStruct.galMask);

fgs=(tCoolStruct.(['in' gasField]).gasMass(tCoolStruct.galMask)+...
    tCoolStruct.(['in' gasField]).sfrMass(tCoolStruct.galMask))./galMass;
fg=fgs./(1+fgs);

%% black hole stuff
bhQM=bhStruct.inGal.cumEngQM(tCoolStruct.galMask);
bhRM=bhStruct.inGal.cumEngRM(tCoolStruct.galMask);
bhMass=bhStruct.inGal.bhMassMax(tCoolStruct.galMask);
%subs.SubhaloMassInRadType(illustris.partTypeNum('bh')+1,tCoolStruct.galMask).*illUnits.massUnit;

%% adress splashback - within the last 5 Gyr at a distance of at least 0.1 R200
if snap==99
    spbMask = mk_splashback_mask('time',5,'both',0.1);
    
else
    spbMask = false(size(tCoolStruct.galMask));
end
%
%
% spb1=splashback.isSat(:,tCoolStruct.galMask);
% spb2=splashback.isFar(:,tCoolStruct.galMask);
% spbMask=(spb1(1,:)==1 & spb2(1,:)==1 ) | ...
%     (spb1(2,:) & spb2(2,:) )|  ...
%     (spb1(3,:)==1 & spb2(3,:)==1 );

galaxyMask=centralMask & gasTemp>0 & bhMass>0 & ~spbMask(tCoolStruct.galMask);
galType='centrals';


nameTag=[gasField '_' galType '_snp' num2str(snap) '_' simDisplayName];




%% get tvir for host fofs
tvir=fofs.Group_T_Crit200(subsInfo.hostFof+1);
tvir=tvir(tCoolStruct.galMask);
global tvirMean
tvirMean=mk_meanMedian_bin(log10(galMass(galaxyMask)),log10(tvir(galaxyMask)),'nb',20);

%% build stellar to halo connection
fofID=subsInfo.hostFof(tCoolStruct.galMask)+1;
hostMass=fofs.Group_M_Crit200(fofID).*illUnits.massUnit;
hostMed=mk_meanMedian_bin(log10(galMass(galaxyMask)),log10(hostMass(galaxyMask)),'bins',9:0.5:12.5);

%hostGasFraction=fofs.GroupMassType(illustris.partTypeNum('gas')+1,fofID)./fofs.GroupMass(fofID);
%fgHost=mk_meanMedian_bin(log10(galMass(galaxyMask)),hostGasFraction(galaxyMask),'bins',9:0.5:12.5);


%% plotting stuff

switch gasField
    case 'Gal'
        gasTag='ISM';
    case {'CGM'}
        gasTag='CGM_{in}';
    case {'Out'}
        gasTag='CGM_{out}';
    case 'Sub'
        gasTag='Subhalo';
end


mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
%tcLab=sprintf('$\\log t_\\mathrm{cool,%s}\\,[\\mathrm{Gyr}]$',gasTag);
tcLab='$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$';
ssfrLab='$\log \mathrm{sSFR}\,[\mathrm{yr^{-1}}]$';
%sfreLab=sprintf('$\\log \\mathrm{SFRe}\\,[\mathrm{yr^{-1}}] $';
% mgasLab=sprintf('$\\log M_\\mathrm{gas,%s}\\,[\\mathrm{M_\\odot}]$',gasTag);
% tempLab=sprintf('$\\log T_{\\mathrm{%s}} \\,[\\mathrm{K}]$',gasTag);
% temp2Lab=sprintf('$\\log T_{\\mathrm{%s}}/T_{\\mathrm{vir}} $',gasTag);
% entLab=sprintf('$\\log S_{\\mathrm{%s}} \\,[\\mathrm{KeV\\,cm^2}]$',gasTag);
% fgsLab=sprintf('$\\log M_\\mathrm{gas,%s}/M_\\mathrm{star}$',gasTag);
% fgLab=sprintf('$\\log M_\\mathrm{gas,%s}/M_\\mathrm{Baryon}$',gasTag);
% %densLab=sprintf('$\\log n_\\mathrm{%s}\\,[\\mathrm{cm^{-3}}]$',gasTag);
% densLab=sprintf('$\\log n_\\mathrm{%s}\\,[\\mathrm{cm^{-3}}]$',gasTag);

mgasLab='$\log M_\mathrm{gas}\,[\mathrm{M_\odot}]$';
tempLab='$\log T \,[\mathrm{K}]$';
temp2Lab='$\log T/T_{\mathrm{vir}} $';
entLab='$\log K \,[\mathrm{keV\,cm^2}]$';
fgsLab='$\log M_\mathrm{gas}/M_\mathrm{star}$';
fgLab='$\log M_\mathrm{gas}/M_\mathrm{Baryons}$';
densLab='$\log n\,[\mathrm{cm^{-3}}]$';

QMLab='$\log E_\mathrm{Thermal}$';
RMLab='$\log E_\mathrm{Kinetic}$';
bhLab='$\log E_\mathrm{Kinetic}/E_\mathrm{Thermal}$';
bhLab2='$\log M_{\mathrm{BH}}\,[\mathrm{M_\odot}]$';
% edissLab=sprintf('$\\log \\dot{E}_\\mathrm{dis,%s}$',gasTag);
% sigmaLab=sprintf('$\\log \\sigma_{\\mathrm{%s}}\\,[\\mathrm{km/sec}]$',gasTag);
% machLab=sprintf('$\\log \\mathcal{M}_\\mathrm{%s}$',gasTag);
% machLab2=sprintf('$\\log \\mathcal{M}_\\mathrm{EW,%s}$',gasTag);

%cmap=brewermap(256,'RdYlBu');
cmap=brewermap(256,'*Spectral');

cmapTemp=cmap;
cmapDens=plasma(256); %   brewermap(256,'*YlGnBu');
%cmapTc=brewermap(256,'YlOrRd');
cmapTc= flipud(viridis(256)); %plasma(256);% magma(256);
cmapBH=cmap; %brewermap(256,'PuBuGn');
cmapSF=flipud(cmapTemp);
cmapRat=plasma(256); %          brewermap(256,'YlGn');
fprintf('ready to plot \n')
%pause

%% Mass vs Temp
if tempFlag
    % create contours
    filt=fspecial('disk',6);
    xdata=log10(galMass(galaxyMask));
    ydata=log10(gasTemp(galaxyMask));
    %obal popCont
    popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
    yl=[4 8];
    xl=[9 12.5];
    
    % create tree
    minLev=7;
    splitParam=30;
    tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
    
    
    % ssfr
    plot_massTemp('tree',xdata,ydata,log10(ssfr(galaxyMask)),tre,ssfrLab,flipud(cmap),xl,yl,popCont,gasField,[-15 -9])
    fname=sprintf('massTemp_%s_%s','ssfr',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % sfre
    %plot_massTemp('tree',xdata,ydata,log10(sfre(galaxyMask)),sfreLab,flipud(cmap))
    
    %fname=sprintf('massTemp_%s_%s','sfre',nameTag);
    %if  printFlag; if  printFlag; printout_fig(gcf,fname,'fig','v'); end end
    
    % tcool
    plot_massTemp('tree',xdata,ydata,log10(tc(galaxyMask)),tre,tcLab,cmapTc,xl,yl,popCont,gasField,[-2.5 1.5])
    
    fname=sprintf('massTemp_%s_%s','tcool',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % gasMass
    plot_massTemp('tree',xdata,ydata,log10(fg(galaxyMask)),tre,fgLab,cmapDens,xl,yl,popCont,gasField,[-3.5 0])
    
    fname=sprintf('massTemp_%s_%s','mf',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % gas density
    plot_massTemp('tree',xdata,ydata,log10(gasDens(galaxyMask)),tre,densLab,cmapDens,xl,yl,popCont,gasField,[-3 -0.5])
    fname=sprintf('massTemp_%s_%s','dens',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % QM
    plot_massTemp('tree',xdata,ydata,log10(bhQM(galaxyMask)),tre,QMLab,cmap,xl,yl,popCont,gasField,[5 9])
    fname=sprintf('massTemp_%s_%s','QM',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % KM
    plot_massTemp('tree',xdata,ydata,log10(bhRM(galaxyMask)),tre,RMLab,cmap,xl,yl,popCont,gasField,[5 9])
    
    fname=sprintf('massTemp_%s_%s','KM',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    %RM/QM
    plot_massTemp('tree',xdata,ydata,log10(bhRM(galaxyMask)./bhQM(galaxyMask)),tre,bhLab,cmap,xl,yl,popCont,gasField,[-5 0])
    
    fname=sprintf('massTemp_%s_%s','bhRat',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    %BH Mass
    plot_massTemp('tree',xdata,ydata,log10(bhMass(galaxyMask)),tre,bhLab2,cmap,xl,yl,popCont,gasField,[6 10])
    fname=sprintf('massTemp_%s_%s','bhMass',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % % Mach
    % plot_massTemp('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap,gasField)
    %
    % fname=sprintf('massTemp_%s_%s','Mach',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % Mach ew
    % plot_massTemp('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap,gasField)
    %
    % fname=sprintf('massTemp_%s_%s','MachEW',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % e dissipation
    % plot_massTemp('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,xl,yl,popCont,gasField,[0 6])
    
    % fname=sprintf('massTemp_%s_%s','eDisp',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % sigma - vel dispersion
    % plot_massTemp('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,xl,yl,popCont,gasField,[1 2.5])
    %
    % fname=sprintf('massTemp_%s_%s','sigma',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
end


%% mass vs entropy

if entFlag
    
    % create contours
    filt=fspecial('disk',6);
    xdata=log10(galMass(galaxyMask));
    ydata=log10(gasEnt(galaxyMask));
    
    popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
    yl=[-3 3];
    xl=[9 12.5];
    
    
    % create quenched contour
    qfilt=fspecial('disk',12);
    xqdata=log10(galMass(galaxyMask & qMask));
    yqdata=log10(gasEnt(galaxyMask & qMask));
    
    qCont=plot_population_contour(xqdata,yqdata,'smooth',qfilt,'noplot');
    
    % create sf contour
    xsfdata=log10(galMass(galaxyMask & ~qMask));
    ysfdata=log10(gasEnt(galaxyMask & ~qMask));
    
    sfCont=plot_population_contour(xsfdata,ysfdata,'smooth',qfilt,'noplot');
    
    
    % create tree
    minLev=7;
    splitParam=30;
    tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
    
    
    
    % ssfr
    plot_massEnt3('tree',xdata,ydata,log10(ssfr(galaxyMask)),tre,ssfrLab,cmapSF,xl,yl,popCont,qCont,sfCont,gasField,[-15 -9],hostMed)
    
    fname=sprintf('massEnt_%s_%s','ssfr',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % sfre
    %plot_massEnt3('tree',xdata,ydata,log10(sfre(galaxyMask)),sfreLab,flipud(cmap))
    %fname=sprintf('massEnt_%s_%s','sfre',nameTag);
    %if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    
    % tcool
    plot_massEnt3('tree',xdata,ydata,log10(tc(galaxyMask)),tre,{tcLab 'in'},cmapTc,xl,yl,popCont,qCont,sfCont,gasField,[-2.5 1.5],hostMed)
    fname=sprintf('massEnt_%s_%s','tcool',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % gasMass
    
    % plot_massEnt3('tree',xdata,ydata,log10(fgs(galaxyMask)),tre,fgsLab,cmap,gasField)
    plot_massEnt3('tree',xdata,ydata,log10(fg(galaxyMask)),tre,{fgLab 'in'},cmapDens,xl,yl,popCont,qCont,sfCont,gasField,[-3.5 0],hostMed)
    fname=sprintf('massEnt_%s_%s','mf',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % gas density
    plot_massEnt3('tree',xdata,ydata,log10(gasDens(galaxyMask)),tre,{densLab 'in'},cmapDens,xl,yl,popCont,qCont,sfCont,gasField,[-4 -1],hostMed)
    fname=sprintf('massEnt_%s_%s','dens',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    %
    %     % QM
    %     plot_massEnt3('tree',xdata,ydata,log10(bhQM(galaxyMask)),tre,QMLab,cmap,xl,yl,popCont,qCont,sfCont,gasField,[5 9],hostMed)
    %     fname=sprintf('massEnt_%s_%s','QM',nameTag);
    %     if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %     % KM
    %     plot_massEnt3('tree',xdata,ydata,log10(bhRM(galaxyMask)),tre,RMLab,cmap,xl,yl,popCont,qCont,sfCont,gasField,[5 9],hostMed)
    %     fname=sprintf('massEnt_%s_%s','KM',nameTag);
    %     if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    
    %RM/QM
    plot_massEnt3('tree',xdata,ydata,log10(bhRM(galaxyMask)./bhQM(galaxyMask)),tre,bhLab,cmapRat,xl,yl,popCont,qCont,sfCont,gasField,[-5 0],hostMed)
    fname=sprintf('massEnt_%s_%s','bhRat',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    %BH Mass
    plot_massEnt3('tree',xdata,ydata,log10(bhMass(galaxyMask)),tre,bhLab2,cmapBH,xl,yl,popCont,qCont,sfCont,gasField,[6 10],hostMed)
    fname=sprintf('massEnt_%s_%s','bhMass',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    % Temperature
    tf=gasTemp(galaxyMask)./tvir(galaxyMask);
    plot_massEnt3('tree',xdata,ydata,log10(tf),tre,{temp2Lab 'in'},cmap,xl,yl,popCont,qCont,sfCont,gasField,[-1.5 1.5],hostMed)
    fname=sprintf('massEnt_%s_%s','temp',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    
    
    
    % % Mach
    % plot_massEnt3('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap,gasField)
    %
    % fname=sprintf('massEnt_%s_%s','Mach',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % Mach ew
    % plot_massEnt3('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap,gasField)
    %
    % fname=sprintf('massEnt_%s_%s','MachEW',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %
    %
    % % e dissipation
    % plot_massEnt3('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,xl,yl,popCont,gasField,[0 6],hostMed)
    %
    % fname=sprintf('massEnt_%s_%s','eDisp',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % sigma - vel dispersion
    % plot_massEnt3('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,xl,yl,popCont,gasField,[1 2.5],hostMed)
    %
    % fname=sprintf('massEnt_%s_%s','sigma',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
end

%% Mass vs. SSFR
if ssfrFlag
    filt=fspecial('disk',6);
    xdata=log10(galMass(galaxyMask));
    ydata=log10(ssfr(galaxyMask));
    popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
    
    minLev=6;
    splitParam=30;
    
    
    xl=[9 12.5];
    yl=[-16.5 -8];
    tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
    
    
    % tcool
    plot_massSsfr2('tree',xdata,ydata,log10(tc(galaxyMask)),tre,tcLab,cmap,popCont,hostMed,[-2.5 1.5],log10(ssfrThresh));
    fname=sprintf('massSsfr_%s_%s','tcool',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    
    % Entropy
    plot_massSsfr2('tree',xdata,ydata,log10(gasEnt(galaxyMask)),tre,entLab,cmap,popCont,hostMed,[-3 3],log10(ssfrThresh));
    fname=sprintf('massSsfr_%s_%s','ent',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    
    % % gasMass
    %
    % plot_massSsfr('tree',xdata,ydata,log10(fgs(galaxyMask)),tre,fgsLab,cmap,[-3 0.5]);
    %
    % fname=sprintf('massSsfr_%s_%s','mf',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % QM
    % plot_massSsfr('tree',xdata,ydata,log10(bhQM(galaxyMask)),tre,QMLab,cmap,[5 9]);
    % fname=sprintf('massSsfr_%s_%s','QM',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % KM
    % plot_massSsfr('tree',xdata,ydata,log10(bhRM(galaxyMask)),tre,RMLab,cmap,[5 9]);
    % fname=sprintf('massSsfr_%s_%s','KM',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %
    %     %RM/QM
    %     plot_massSsfr('tree',xdata,ydata,log10(bhRM(galaxyMask)./bhQM(galaxyMask)),tre,bhLab,cmap,[-5 0]);
    %     fname=sprintf('massSsfr_%s_%s','bhRat',nameTag);
    %     if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %
    %     %BH Mass
    %     plot_massSsfr('tree',xdata,ydata,log10(bhMass),tre,bhLab2,cmap,[6 10])
    %     fname=sprintf('massSsfr_%s_%s','bhMass',nameTag);
    %     if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %
    % % Mach
    % plot_massSsfr('tree',xdata,ydata,log10(gasMach(galaxyMask)),tre,machLab,cmap);
    %
    % fname=sprintf('massSsfr_%s_%s','Mach',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % Mach ew
    % plot_massSsfr('tree',xdata,ydata,log10(gasMach2(galaxyMask)),tre,machLab2,cmap)
    %
    % fname=sprintf('massSsfr_%s_%s','MachEW',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %
    %
    % % e dissipation
    % plot_massSsfr('tree',xdata,ydata,log10(gasEDisp(2,galaxyMask)),tre,edissLab,cmap,[0 6]);
    %
    % fname=sprintf('massSsfr_%s_%s','eDisp',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    % % sigma - vel dispersion
    % plot_massSsfr('tree',xdata,ydata,log10(gasVDisp(1,galaxyMask)),tre,sigmaLab,cmap,[1 2.5]);
    %
    % fname=sprintf('massSsfr_%s_%s','sigma',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %
    % % temp
    % plot_massSsfr('tree',xdata,ydata,log10(gasTemp(galaxyMask)),tre,tempLab,cmap,[4 8]);
    %
    % fname=sprintf('massSsfr_%s_%s','Temp',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %
    % % Ent
    % plot_massSsfr('tree',xdata,ydata,log10(gasEnt(galaxyMask)),tre,entLab,cmap,[-3 2.5]);
    %
    % fname=sprintf('massSsfr_%s_%s','Ent',nameTag);
    % if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %
    %
    %
    %
    %
end

%% bhrat vs. entropy
if bhratFlag
    
    base=10^0.5*1e-8;
    scat=0.5;
    bhrat=bhRM(galaxyMask)./bhQM(galaxyMask);
    msk=bhrat==0 | isinf(bhrat);
    bhrat(msk)=base.*10.^(scat.*rand(1,sum(msk)));
    
    
    filt=fspecial('disk',6);
    xdata=log10(bhrat);
    ydata=log10(gasEnt(galaxyMask));
    
    popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
    yl=[-3 3];
    xl=[-7.6 0];
    
    % create quenched contour
    qfilt=fspecial('disk',12);
    brq=bhRM(galaxyMask & qMask)./bhQM(galaxyMask & qMask);
    msk=brq==0 | isinf(brq);
    brq(msk)=base.*10.^(scat.*rand(1,sum(msk)));
    
    
    xqdata=log10(brq);
    yqdata=log10(gasEnt(galaxyMask & qMask));
    
    qCont=plot_population_contour(xqdata,yqdata,'smooth',qfilt,'noplot');
    
    % create sf contour
    brs=bhRM(galaxyMask & ~qMask)./bhQM(galaxyMask & ~qMask);
    msk=brs==0 | isinf(brs);
    brs(msk)=base.*10.^(scat.*rand(1,sum(msk)));
    
    xsfdata=log10(brs);
    ysfdata=log10(gasEnt(galaxyMask & ~qMask));
    
    sfCont=plot_population_contour(xsfdata,ysfdata,'smooth',qfilt,'noplot');
    
    
    minLev=7;
    splitParam=30;
    
    tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);
    
    % ssfr
    plot_bhratEnt2('tree',xdata,ydata,log10(ssfr(galaxyMask)),tre,ssfrLab,cmapSF,xl,yl,popCont,qCont,sfCont,'Gal',[-16.5 -9],hostMed);
    
    fname=sprintf('bhratEnt_%s_%s','ssfr',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    
end
