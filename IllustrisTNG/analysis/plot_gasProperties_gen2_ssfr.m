%% Plot results for gas properties in Centrals in TNG
%% load data
snap=99;
global simDisplayName
if readFlag
    
    global DEFAULT_MATFILE_DIR
    %load([DEFAULT_MATFILE_DIR '/cooling_times_z0_' simDisplayName '.mat'])
    load([DEFAULT_MATFILE_DIR '/gasProperties_snp' num2str(snap) '_' simDisplayName '.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_snp' num2str(snap) '_' simDisplayName '.mat'])
    
    
    %     bhStruct.inGal=illustris.utils.read_catalog('bh_inGal','folder','bhProps');
    %     tCoolStruct.inGal=illustris.utils.read_catalog('gasProps_inGal','folder','gasProperties');
    
    
    %     global DRACOFLAG
    %     if DRACOFLAG
    %         fofs=illustris.groupcat.loadHalos(bp,99);
    %         subs=illustris.groupcat.loadSubhalos(bp,99);
    %
    %     else
    %         load([DEFAULT_MATFILE_DIR '/tng100_z0_fofs_subs.mat'])
    %     end
    %
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
    
end
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
fofs=illustris.utils.addTvirFofs(fofs);
illustris.utils.set_illUnits(snap)
global illUnits

%% define what we are plotting
% identify centrals
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

%gasField='Gal';

%global simDisplayName

% get usefel stuff
galMass=tCoolStruct.galMass(tCoolStruct.galMask);  % galaxy stellar mass

galGasMass=tCoolStruct.inGal.gasMass(:,tCoolStruct.galMask);
cgmGasMass=tCoolStruct.inCGM.gasMass(:,tCoolStruct.galMask);
outGasMass=tCoolStruct.inOut.gasMass(:,tCoolStruct.galMask);

sfr=subs.SubhaloSFRinRad(tCoolStruct.galMask);  % sfr in galaxy
sfrBase=10.^0.75.*1e-15;
ssfr=illustris.utils.calc_ssfr(subs,'base',sfrBase,'scatter',0.25);% sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));
ssfr(ssfr<sfrBase)=sfrBase;

%ssfr=illustris.utils.calc_ssfr(subs);

ssfr=ssfr(tCoolStruct.galMask);

ssfrThresh=1e-11;
qMask=ssfr<=ssfrThresh;

%% gas properties
galTc=tCoolStruct.inGal.meanTcMW(1,tCoolStruct.galMask);
cgmTc=tCoolStruct.inCGM.meanTcMW(1,tCoolStruct.galMask);
outTc=tCoolStruct.inOut.meanTcMW(1,tCoolStruct.galMask);

galTemp=tCoolStruct.inGal.meanTempMW(1,tCoolStruct.galMask);
cgmTemp=tCoolStruct.inCGM.meanTempMW(1,tCoolStruct.galMask);
outTemp=tCoolStruct.inOut.meanTempMW(1,tCoolStruct.galMask);

galEnt=tCoolStruct.inGal.meanEntMW(1,tCoolStruct.galMask);
cgmEnt=tCoolStruct.inCGM.meanEntMW(1,tCoolStruct.galMask);
outEnt=tCoolStruct.inOut.meanEntMW(1,tCoolStruct.galMask);

galDens=tCoolStruct.inGal.meanDensN(1,tCoolStruct.galMask);
cgmDens=tCoolStruct.inCGM.meanDensN(1,tCoolStruct.galMask);
outDens=tCoolStruct.inOut.meanDensN(1,tCoolStruct.galMask);

fgs=(tCoolStruct.inGal.gasMass(tCoolStruct.galMask)+...
    tCoolStruct.inGal.sfrMass(tCoolStruct.galMask))./galMass;
fg=fgs./(1+fgs);

%% black hole stuff
bhQM=bhStruct.inGal.cumEngQM(tCoolStruct.galMask);
bhRM=bhStruct.inGal.cumEngRM(tCoolStruct.galMask);
bhMass=bhStruct.inGal.bhMassMax(tCoolStruct.galMask);

spbMask = mk_splashback_mask('time',5,'both',0.1);


galaxyMaskBase=centralMask & galTemp>0 & bhMass>0 & ~spbMask(tCoolStruct.galMask);
galType='centrals';





%% get tvir for host fofs
tvir=fofs.Group_T_Crit200(subsInfo.hostFof+1);
tvir=tvir(tCoolStruct.galMask);
% global tvirMean
% tvirMean=mk_meanMedian_bin(log10(galMass(galaxyMask)),log10(tvir(galaxyMask)),'nb',20);

%% build stellar to halo connection
hostMass=fofs.Group_M_Crit200(subsInfo.hostFof(tCoolStruct.galMask)+1).*illUnits.massUnit;
%hostMed=mk_meanMedian_bin(log10(galMass(galaxyMask)),log10(hostMass(galaxyMask)),'bins',9:0.5:12.5);

%% plotting stuff

% switch gasField
%     case 'Gal'
%         ismTag='ISM';
%     case {'CGM'}
%         gasTag='CGM_{in}';
%         case {'Out'}
%         gasTag='CGM_{out}';
%     case 'Sub'
%         gasTag='Subhalo';
% end
%

mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
%tcLab=sprintf('$\\log t_\\mathrm{cool,%s}\\,[\\mathrm{Gyr}]$',gasTag);
tcLab='$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$';
ssfrLab='$\log \mathrm{sSFR}\,[\mathrm{yr^{-1}}]$';
%sfreLab=sprintf('$\\log \\mathrm{SFRe}\\,[\mathrm{yr^{-1}}] $';
%mgasLab=sprintf('$\\log M_\\mathrm{gas,%s}\\,[\\mathrm{M_\\odot}]$',gasTag);
tempLab=sprintf('$\\log T\\,[\\mathrm{K}]$');
temp2Lab=sprintf('$\\log T/T_{\\mathrm{vir}} $');
entLab=sprintf('$ S \\,[\\mathrm{KeV\\,cm^2}]$');
%fgsLab=sprintf('$\\log M_\\mathrm{gas,%s}/M_\\mathrm{star}$',gasTag);
fgLab=sprintf('$\\log M_\\mathrm{gas,%s}/M_\\mathrm{Baryon}$','ISM');
%densLab=sprintf('$\\log n_\\mathrm{%s}\\,[\\mathrm{cm^{-3}}]$',gasTag);
densLab=sprintf('$\\log n\\,[\\mathrm{cm^{-3}}]$');
QMLab='$\log E_\mathrm{HAM}$';
RMLab='$\log E_\mathrm{LAM}$';
bhLab='$\log E_\mathrm{LAM}/E_\mathrm{HAM}$';
bhLab2='$\log M_{\mathrm{BH}}\,[\mathrm{M_\odot}]$';
% edissLab=sprintf('$\\log \\dot{E}_\\mathrm{dis,%s}$',gasTag);
% sigmaLab=sprintf('$\\log \\sigma_{\\mathrm{%s}}\\,[\\mathrm{km/sec}]$',gasTag);
% machLab=sprintf('$\\log \\mathcal{M}_\\mathrm{%s}$',gasTag);
% machLab2=sprintf('$\\log \\mathcal{M}_\\mathrm{EW,%s}$',gasTag);

%cmap=brewermap(256,'RdYlBu');
cmap=brewermap(256,'*Spectral');
cc=brewermap(8,'Set1');
cmapTemp=cmap;
cmapDens=brewermap(256,'*YlGnBu');
cmapTc=brewermap(256,'YlOrRd');
cmapBH=brewermap(256,'PuBuGn');
cmapSF=flipud(cmapTemp);
cmapRat=brewermap(256,'YlGn');
fprintf('ready to plot \n')
%pause


for i=5 % 1:5
    switch(i)
        case 1
            massMask=galMass >= 1e9 & galMass <1e10;
            mtag='$[10^9,10^{10}]$';
            ptag='m9';
        case 2
            massMask=galMass >= 1e10 & galMass <1e11;
            mtag='$[10^{10},10^{11}]$';
            ptag='m10';
        case 3
            massMask=galMass >= 1e11 ;
            mtag='$\ge10^{11}$';
            ptag='m11';
        case 4
            massMask=galMass >= 1e9 ;
            mtag='$\ge10^9$';
            ptag='mAll';
            
        case 5
            massMask=log10(galMass)>=10.3 & log10(galMass)<10.7;
            mtag='$[10^{10.3},10^{10.7}]$';
            ptag='m1037';
    end
    
    galaxyMask=galaxyMaskBase & massMask;
    
    nameTag=['comp3_' ptag '_' galType '_snp' num2str(snap) '_' simDisplayName];
    
    %% plot vs. ssfr
    
    xdata=log10(ssfr(galaxyMask));
    bins=-14.5:0.75:-8.5;
    xl=[-14.5 -8.5];
    
    %% tcool
    yl=[-3 2];
    yGal=mk_meanMedian_bin(xdata,galTc(galaxyMask),'bins',bins);
    yCGM=mk_meanMedian_bin(xdata,cgmTc(galaxyMask),'bins',bins);
    yOut=mk_meanMedian_bin(xdata,outTc(galaxyMask),'bins',bins);
    
    mmG=yGal.binCount>0;
    galPolY=log10([yGal.yQuarts(2,mmG) fliplr(yGal.yQuarts(3,mmG))]);
    galPolX=[yGal.xMedian(mmG) fliplr(yGal.xMedian(mmG))];
    galPolY(isinf(galPolY))=-10;
    
    mmC=yCGM.binCount>0;
    cgmPolY=log10([yCGM.yQuarts(2,mmC) fliplr(yCGM.yQuarts(3,mmC))]);
    cgmPolX=[yCGM.xMedian(mmC) fliplr(yCGM.xMedian(mmC))];
    cgmPolY(isinf(cgmPolY))=-10;
    
    mmO=yOut.binCount>0;
    outPolY=log10([yOut.yQuarts(2,mmO) fliplr(yOut.yQuarts(3,mmO))]);
    outPolX=[yOut.xMedian(mmO) fliplr(yOut.xMedian(mmO))];
    outPolY(isinf(outPolY))=-10;
    
    figure('position',[1432 421 1000 750],'Color','w')
    yyaxis right
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.k','linewidth',1)
    ylabelmine('$N/N_\mathrm{Tot}$');;
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian(mmG),log10(yGal.yMedian(mmG)),'-','color',cc(1,:),'DisplayName','ISM','linewidth',2);
    hold on
    h(2)=plot(yCGM.xMedian(mmC),log10(yCGM.yMedian(mmC)),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',2);
    h(3)=plot(yOut.xMedian(mmO),log10(yOut.yMedian(mmO)),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',2);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
    
    nTag=sprintf('%s, %i',mtag,sum(galaxyMask));
    xfac=0.7; yfac=1.02; 
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',20,'fontweight','bold')
    grid
    hl=legend(h);set(hl,'Interpreter','latex','fontsize',16,'location','SouthWest');
    
    xlim(xl);ylim(yl);
    xlabelmine(ssfrLab);
    ylabelmine(tcLab);
    set(gca,'Fontsize',18)
    
    fname=sprintf('ssfrTcool_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    
    %% temp
    yl=[4 7];
    yGal=mk_meanMedian_bin(xdata,galTemp(galaxyMask),'bins',bins);
    yCGM=mk_meanMedian_bin(xdata,cgmTemp(galaxyMask),'bins',bins);
    yOut=mk_meanMedian_bin(xdata,outTemp(galaxyMask),'bins',bins);
    
    mmG=yGal.binCount>0;
    galPolY=log10([yGal.yQuarts(2,mmG) fliplr(yGal.yQuarts(3,mmG))]);
    galPolX=[yGal.xMedian(mmG) fliplr(yGal.xMedian(mmG))];
    galPolY(isinf(galPolY))=-10;
    
    mmC=yCGM.binCount>0;
    cgmPolY=log10([yCGM.yQuarts(2,mmC) fliplr(yCGM.yQuarts(3,mmC))]);
    cgmPolX=[yCGM.xMedian(mmC) fliplr(yCGM.xMedian(mmC))];
    cgmPolY(isinf(cgmPolY))=-10;
    
    mmO=yOut.binCount>0;
    outPolY=log10([yOut.yQuarts(2,mmO) fliplr(yOut.yQuarts(3,mmO))]);
    outPolX=[yOut.xMedian(mmO) fliplr(yOut.xMedian(mmO))];
    outPolY(isinf(outPolY))=-10;
    
    figure('position',[1432 421 1000 750],'Color','w')
    yyaxis right
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.k','linewidth',1)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian,log10(yGal.yMedian),'-','color',cc(1,:),'DisplayName','ISM','linewidth',2);
    hold on
    h(2)=plot(yCGM.xMedian,log10(yCGM.yMedian),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',2);
    h(3)=plot(yOut.xMedian,log10(yOut.yMedian),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',2);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
    
    nTag=sprintf('%s, %i',mtag,sum(galaxyMask));
    xfac=0.7; yfac=1.02; 
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',20,'fontweight','bold')
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',16,'location','SouthWest');
    
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,18);
    ylabelmine(tempLab,18);
    
    set(gca,'Fontsize',18)
    
    fname=sprintf('ssfrTemp_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    %% temp / tvir
    yl=[-1.5 0.5];
    yGal=mk_meanMedian_bin(xdata,galTemp(galaxyMask)./tvir(galaxyMask),'bins',bins);
    yCGM=mk_meanMedian_bin(xdata,cgmTemp(galaxyMask)./tvir(galaxyMask),'bins',bins);
    yOut=mk_meanMedian_bin(xdata,outTemp(galaxyMask)./tvir(galaxyMask),'bins',bins);
    
    mmG=yGal.binCount>0;
    galPolY=log10([yGal.yQuarts(2,mmG) fliplr(yGal.yQuarts(3,mmG))]);
    galPolX=[yGal.xMedian(mmG) fliplr(yGal.xMedian(mmG))];
    galPolY(isinf(galPolY))=-10;
    
    mmC=yCGM.binCount>0;
    cgmPolY=log10([yCGM.yQuarts(2,mmC) fliplr(yCGM.yQuarts(3,mmC))]);
    cgmPolX=[yCGM.xMedian(mmC) fliplr(yCGM.xMedian(mmC))];
    cgmPolY(isinf(cgmPolY))=-10;
    
    mmO=yOut.binCount>0;
    outPolY=log10([yOut.yQuarts(2,mmO) fliplr(yOut.yQuarts(3,mmO))]);
    outPolX=[yOut.xMedian(mmO) fliplr(yOut.xMedian(mmO))];
    outPolY(isinf(outPolY))=-10;
    
    figure('position',[1432 421 1000 750],'Color','w')
    yyaxis right
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.k','linewidth',1)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian,log10(yGal.yMedian),'-','color',cc(1,:),'DisplayName','ISM','linewidth',2);
    hold on
    h(2)=plot(yCGM.xMedian,log10(yCGM.yMedian),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',2);
    h(3)=plot(yOut.xMedian,log10(yOut.yMedian),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',2);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
    nTag=sprintf('%s, %i',mtag,sum(galaxyMask));
    xfac=0.7; yfac=1.02; 
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',20,'fontweight','bold')
   
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',16,'location','SouthWest');
    
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,18);
    ylabelmine(temp2Lab,18);
    
    set(gca,'Fontsize',18)
    
    fname=sprintf('ssfrTemp2_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    %% density
    yl=[-4.5 -0.5];
    yGal=mk_meanMedian_bin(xdata,galDens(galaxyMask),'bins',bins);
    yCGM=mk_meanMedian_bin(xdata,cgmDens(galaxyMask),'bins',bins);
    yOut=mk_meanMedian_bin(xdata,outDens(galaxyMask),'bins',bins);
    
     mmG=yGal.binCount>0;
    galPolY=log10([yGal.yQuarts(2,mmG) fliplr(yGal.yQuarts(3,mmG))]);
    galPolX=[yGal.xMedian(mmG) fliplr(yGal.xMedian(mmG))];
    galPolY(isinf(galPolY))=-10;
    
    mmC=yCGM.binCount>0;
    cgmPolY=log10([yCGM.yQuarts(2,mmC) fliplr(yCGM.yQuarts(3,mmC))]);
    cgmPolX=[yCGM.xMedian(mmC) fliplr(yCGM.xMedian(mmC))];
    cgmPolY(isinf(cgmPolY))=-10;
    
    mmO=yOut.binCount>0;
    outPolY=log10([yOut.yQuarts(2,mmO) fliplr(yOut.yQuarts(3,mmO))]);
    outPolX=[yOut.xMedian(mmO) fliplr(yOut.xMedian(mmO))];
    outPolY(isinf(outPolY))=-10;
    
    
    figure('position',[1432 421 1000 750],'Color','w')
    yyaxis right
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.k','linewidth',1)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian,log10(yGal.yMedian),'-','color',cc(1,:),'DisplayName','ISM','linewidth',2);
    hold on
    h(2)=plot(yCGM.xMedian,log10(yCGM.yMedian),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',2);
    h(3)=plot(yOut.xMedian,log10(yOut.yMedian),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',2);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
    nTag=sprintf('%s, %i',mtag,sum(galaxyMask));
    xfac=0.7; yfac=1.02; 
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',20,'fontweight','bold')
   
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',16,'location','SouthWest');
    
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,18);
    ylabelmine(densLab,18);
    
    set(gca,'Fontsize',18)
    fname=sprintf('ssfrDens_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    %% gas mass fraction density
    yl=[-3 0];
    yGal=mk_meanMedian_bin(xdata,fg(galaxyMask),'bins',bins);
    
    
    mmG=yGal.binCount>0;
    galPolY=log10([yGal.yQuarts(2,mmG) fliplr(yGal.yQuarts(3,mmG))]);
    galPolX=[yGal.xMedian(mmG) fliplr(yGal.xMedian(mmG))];
    galPolY(isinf(galPolY))=-10;
    
    
    figure('position',[1432 421 1000 750],'Color','w')
    yyaxis right
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.k','linewidth',1)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian,log10(yGal.yMedian),'-','color',cc(1,:),'DisplayName','ISM','linewidth',2);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    nTag=sprintf('%s, %i',mtag,sum(galaxyMask));
    xfac=0.7; yfac=1.02; 
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',20,'fontweight','bold')
   
    grid
    
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,18);
    ylabelmine(fgLab,18);
    
    set(gca,'Fontsize',18)
    fname=sprintf('ssfrgasFrac_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
  
    %% entropy
    yl=[-2 2.5];
    yGal=mk_meanMedian_bin(xdata,galEnt(galaxyMask),'bins',bins);
    yCGM=mk_meanMedian_bin(xdata,cgmEnt(galaxyMask),'bins',bins);
    yOut=mk_meanMedian_bin(xdata,outEnt(galaxyMask),'bins',bins);
    
     mmG=yGal.binCount>0;
    galPolY=log10([yGal.yQuarts(2,mmG) fliplr(yGal.yQuarts(3,mmG))]);
    galPolX=[yGal.xMedian(mmG) fliplr(yGal.xMedian(mmG))];
    galPolY(isinf(galPolY))=-10;
    
    mmC=yCGM.binCount>0;
    cgmPolY=log10([yCGM.yQuarts(2,mmC) fliplr(yCGM.yQuarts(3,mmC))]);
    cgmPolX=[yCGM.xMedian(mmC) fliplr(yCGM.xMedian(mmC))];
    cgmPolY(isinf(cgmPolY))=-10;
    
    mmO=yOut.binCount>0;
    outPolY=log10([yOut.yQuarts(2,mmO) fliplr(yOut.yQuarts(3,mmO))]);
    outPolX=[yOut.xMedian(mmO) fliplr(yOut.xMedian(mmO))];
    outPolY(isinf(outPolY))=-10;
    
    
    figure('position',[1432 421 1000 750],'Color','w')
    
    yyaxis right
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.k','linewidth',1)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian,log10(yGal.yMedian),'-','color',cc(1,:),'DisplayName','ISM','linewidth',2);
    hold on
    h(2)=plot(yCGM.xMedian,log10(yCGM.yMedian),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',2);
    h(3)=plot(yOut.xMedian,log10(yOut.yMedian),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',2);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
    nTag=sprintf('%s, %i',mtag,sum(galaxyMask));
    xfac=0.7; yfac=1.02; 
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',20,'fontweight','bold')
   
    grid
    hl=legend(h);set(hl,'Interpreter','latex','fontsize',16,'location','SouthWest');
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,18);
    ylabelmine(entLab,18);
    set(gca,'Fontsize',18)
    
    fname=sprintf('ssfrEnt_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
    
    %close all
end