%% Plot results for gas properties in Centrals in TNG
%% load data
snap=99;
global DRACOFLAG
global simDisplayName
if readFlag
    
    global DEFAULT_MATFILE_DIR
    %load([DEFAULT_MATFILE_DIR '/cooling_times_z0_' simDisplayName '.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_snp' num2str(snap) '_' simDisplayName '.mat'])
     load([DEFAULT_MATFILE_DIR '/gasProperties_snp' num2str(snap) '_' simDisplayName '.mat'])
 
    
    %     bhStruct.inGal=illustris.utils.read_catalog('bh_inGal','folder','bhProps');
    %     tCoolStruct.inGal=illustris.utils.read_catalog('gasProps_inGal','folder','gasProperties');
       
    if DRACOFLAG
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
        
    else
        loadFofSubTNG100
    end
    
    fofs=illustris.utils.addTvirFofs(fofs);
end
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
illustris.utils.set_illUnits(snap)
global illUnits

%% define what we are plotting
% identify centrals
%centralMask= subsInfo.isCentral(tCoolStruct.galMask);
spbMask = mk_splashback_mask('time',5,'both',0.1);

totMask=subsInfo.isCentral & tCoolStruct.galMask & ~spbMask;



%gasField='Gal';

%global simDisplayName

% get usefel stuff
galMass=tCoolStruct.galMass(totMask);  % galaxy stellar mass

galGasMass=tCoolStruct.inGal.gasMass(:,totMask);
cgmGasMass=tCoolStruct.inCGM.gasMass(:,totMask);
outGasMass=tCoolStruct.inOut.gasMass(:,totMask);

sfr=subs.SubhaloSFRinRad(totMask);  % sfr in galaxy
sfrBase=10.^0.75.*1e-15;
ssfr=illustris.utils.calc_ssfr(subs,'base',sfrBase,'scatter',0.25);% sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));
ssfr(ssfr<sfrBase)=sfrBase;

%ssfr=illustris.utils.calc_ssfr(subs);

ssfr=ssfr(totMask);

ssfrThresh=1e-11;
qMask=ssfr<=ssfrThresh;

%% gas properties
galTc=tCoolStruct.inGal.meanTcMW(1,totMask);
cgmTc=tCoolStruct.inCGM.meanTcMW(1,totMask);
outTc=tCoolStruct.inOut.meanTcMW(1,totMask);

galTemp=tCoolStruct.inGal.meanTempMW(1,totMask);
cgmTemp=tCoolStruct.inCGM.meanTempMW(1,totMask);
outTemp=tCoolStruct.inOut.meanTempMW(1,totMask);

galEnt=tCoolStruct.inGal.meanEntMW(1,totMask);
cgmEnt=tCoolStruct.inCGM.meanEntMW(1,totMask);
outEnt=tCoolStruct.inOut.meanEntMW(1,totMask);

galDens=tCoolStruct.inGal.meanDensN(1,totMask);
cgmDens=tCoolStruct.inCGM.meanDensN(1,totMask);
outDens=tCoolStruct.inOut.meanDensN(1,totMask);

fgs=(tCoolStruct.inGal.gasMass(totMask)+...
    tCoolStruct.inGal.sfrMass(totMask))./galMass;
fg=fgs./(1+fgs);

%% black hole stuff
bhQM=bhStruct.inGal.cumEngQM(totMask);
bhRM=bhStruct.inGal.cumEngRM(totMask);
bhMass=bhStruct.inGal.bhMassMax(totMask);

galaxyMaskBase=true(size(galMass));
galType='centrals';





%% get tvir for host fofs
tvir=fofs.Group_T_Crit200(subsInfo.hostFof+1);
tvir=tvir(totMask);
% global tvirMean
% tvirMean=mk_meanMedian_bin(log10(galMass(galaxyMask)),log10(tvir(galaxyMask)),'nb',20);

%% build stellar to halo connection
hostMass=fofs.Group_M_Crit200(subsInfo.hostFof(totMask)+1).*illUnits.massUnit;
%hostMed=mk_meanMedian_bin(log10(galMass(galaxyMask)),log10(hostMass(galaxyMask)),'bins',9:0.5:12.5);

%% plotting stuff

if ~DRACOFLAG
    tagFont=16;
    axFont=18;
    bigTagFont=18;
    lw=2;
else
    tagFont=30;
    axFont=30;
    bigTagFont=34;
    lw=5;
end


sfTag='Star-Forming';
qTag='Quenched';

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
cc2=brewermap(8,'Accent');
cmapTemp=cmap;
cmapDens=brewermap(256,'*YlGnBu');
cmapTc=brewermap(256,'YlOrRd');
cmapBH=brewermap(256,'PuBuGn');
cmapSF=flipud(cmapTemp);
cmapRat=brewermap(256,'YlGn');
fprintf('ready to plot \n')

ccc=[0.7 0.7 0.7];
ccc2=ccc-0.3;

%pause


for i= 2
    switch(i)
        case 1
            massMask=galMass >= 1e9 & galMass <1e10;
            mtag='$[10^9,10^{10}]$';
            ptag='m9';
            subDir='comp3/gasProp_vs_ssfr/massRange_9_10';
        case 2
            massMask=galMass >= 1e10 & galMass <1e11;
            mtag='$[10^{10},10^{11}]$';
            ptag='m10';
            subDir='comp3/gasProp_vs_ssfr/massRange_10_11';
        case 3
            massMask=galMass >= 1e11 ;
            mtag='$\ge10^{11}$';
            ptag='m11';
            subDir='comp3/gasProp_vs_ssfr/massRange_11_andUp';
        case 4
            massMask=galMass >= 1e9 ;
            mtag='$\ge10^9$';
            ptag='mAll';
            subDir='comp3/gasProp_vs_ssfr/massRange_all';            
        case 5
            massMask=log10(galMass)>=9.5 & log10(galMass)<10.3;
            mtag='$[10^{9.5},10^{10.3}]$';
            ptag='m95';
            subDir='comp3/gasProp_vs_ssfr/massRange_9.5_10.3';
        case 6
            massMask=log10(galMass)>=10.3 & log10(galMass)<10.7;
            mtag='$[10^{10.3},10^{10.7}]$';
            ptag='m1037';
            subDir='comp3/gasProp_vs_ssfr/massRange_10.3_10.7';
        case 7
            massMask=log10(galMass)>=10.7 & log10(galMass)<11.5;
            mtag='$[10^{10.7},10^{11.5}]$';
            ptag='m107';
            subDir='comp3/gasProp_vs_ssfr/massRange_10.7_11.5';
           
            
    end
    
    galaxyMask=galaxyMaskBase & massMask;
    
    nameTag=['comp3_' ptag '_' galType '_snp' num2str(snap) '_' simDisplayName];
    
    %% plot vs. ssfr
    
    xdata=log10(ssfr(galaxyMask));
    bins=linspace(-14.5,-8.5,11);  %-14.5:0.5:-8.5;
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
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.','linewidth',lw,'color',ccc)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    set(gca,'YColor',ccc2)
    
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian(mmG),log10(yGal.yMedian(mmG)),'-','color',cc(1,:),'DisplayName','Galactic Gas','linewidth',lw);
    hold on
    h(2)=plot(yCGM.xMedian(mmC),log10(yCGM.yMedian(mmC)),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',lw);
    h(3)=plot(yOut.xMedian(mmO),log10(yOut.yMedian(mmO)),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',lw);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
    
    nTag=sprintf('%s, %s, %i',simDisplayName,mtag,sum(galaxyMask));
    xfac=0.6; yfac=1.025;
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',tagFont,'fontweight','bold')
   
    % annotations 
    plot(log10(ssfrThresh).*[1 1],yl,'k:','linewidth',lw);
    dx=diff(xl)/30;
    dy=diff(yl)/30;
        
    sfTagPos=[log10(ssfrThresh)+dx yl(2)-dy];
    
    text(sfTagPos(1),sfTagPos(2),sfTag,'color',[0 0 1],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    qTagPos =[log10(ssfrThresh)-6.6*dx yl(2)-dy];
    text(qTagPos(1),qTagPos(2),qTag,'color',[1 0 0],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    
    arrowLen=4*dx;
    sfAstart=[sfTagPos(1) sfTagPos(2)-dy];
    sfAend=[sfTagPos(1)+arrowLen sfTagPos(2)-dy];
    arrow(sfAstart,sfAend,'width',lw,'FaceColor',[0 0 1],'EdgeColor',[ 0 0 1])
    
    qAstart=[log10(ssfrThresh)-dx qTagPos(2)-dy];
    qAend=[qAstart(1)-arrowLen qTagPos(2)-dy];
    arrow(qAstart,qAend,'width',lw,'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
        
    grid
    hl=legend(h);set(hl,'Interpreter','latex','fontsize',tagFont,'location','SouthWest');
    
    xlim(xl);ylim(yl);
    xlabelmine(ssfrLab,tagFont);
    ylabelmine(tcLab,tagFont);
    set(gca,'Fontsize',axFont,'YColor','k')
    
    fname=sprintf('ssfrTcool_%s',nameTag);
    
    if  printFlag; printout_fig(gcf,fname,'subdir',subDir,'fig','v'); end
    
    
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
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.','linewidth',lw,'color',ccc)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    set(gca,'YColor',ccc2)
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian(mmG),log10(yGal.yMedian(mmG)),'-','color',cc(1,:),'DisplayName','Galactic Gas','linewidth',lw);
    hold on
    h(2)=plot(yCGM.xMedian(mmC),log10(yCGM.yMedian(mmC)),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',lw);
    h(3)=plot(yOut.xMedian(mmO),log10(yOut.yMedian(mmO)),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',lw);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
    
    nTag=sprintf('%s, %s, %i',simDisplayName,mtag,sum(galaxyMask));
    xfac=0.6; yfac=1.025;
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',tagFont,'fontweight','bold')
    
    
    % annotations 
     plot(log10(ssfrThresh).*[1 1],yl,'k:','linewidth',lw);
    dx=diff(xl)/30;
    dy=diff(yl)/30;
        
    sfTagPos=[log10(ssfrThresh)+dx yl(2)-dy];
    
    text(sfTagPos(1),sfTagPos(2),sfTag,'color',[0 0 1],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    qTagPos =[log10(ssfrThresh)-6.6*dx yl(2)-dy];
    text(qTagPos(1),qTagPos(2),qTag,'color',[1 0 0],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    
    arrowLen=4*dx;
    sfAstart=[sfTagPos(1) sfTagPos(2)-dy];
    sfAend=[sfTagPos(1)+arrowLen sfTagPos(2)-dy];
    arrow(sfAstart,sfAend,'width',lw,'FaceColor',[0 0 1],'EdgeColor',[ 0 0 1])
    
    qAstart=[log10(ssfrThresh)-dx qTagPos(2)-dy];
    qAend=[qAstart(1)-arrowLen qTagPos(2)-dy];
    arrow(qAstart,qAend,'width',lw,'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
      
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',tagFont,'location','SouthWest');
    
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,tagFont);
    ylabelmine(tempLab,tagFont);
    
        set(gca,'Fontsize',axFont,'YColor','k') 
    
    fname=sprintf('ssfrTemp_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'subdir',subDir,'fig','v'); end
    
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
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.','linewidth',lw,'color',ccc)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    set(gca,'YColor',ccc2)
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian(mmG),log10(yGal.yMedian(mmG)),'-','color',cc(1,:),'DisplayName','Galactic Gas','linewidth',lw);
    hold on
    h(2)=plot(yCGM.xMedian(mmC),log10(yCGM.yMedian(mmC)),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',lw);
    h(3)=plot(yOut.xMedian(mmO),log10(yOut.yMedian(mmO)),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',lw);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
   
     nTag=sprintf('%s, %s, %i',simDisplayName,mtag,sum(galaxyMask));
    xfac=0.6; yfac=1.025;
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',tagFont,'fontweight','bold')
    
    % annotations 
     plot(log10(ssfrThresh).*[1 1],yl,'k:','linewidth',lw);
    dx=diff(xl)/30;
    dy=diff(yl)/30;
        
    sfTagPos=[log10(ssfrThresh)+dx yl(2)-dy];
    
    text(sfTagPos(1),sfTagPos(2),sfTag,'color',[0 0 1],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    qTagPos =[log10(ssfrThresh)-6.6*dx yl(2)-dy];
    text(qTagPos(1),qTagPos(2),qTag,'color',[1 0 0],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    
    arrowLen=4*dx;
    sfAstart=[sfTagPos(1) sfTagPos(2)-dy];
    sfAend=[sfTagPos(1)+arrowLen sfTagPos(2)-dy];
    arrow(sfAstart,sfAend,'width',lw,'FaceColor',[0 0 1],'EdgeColor',[ 0 0 1])
    
    qAstart=[log10(ssfrThresh)-dx qTagPos(2)-dy];
    qAend=[qAstart(1)-arrowLen qTagPos(2)-dy];
    arrow(qAstart,qAend,'width',lw,'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
   
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',tagFont,'location','SouthWest');
    
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,tagFont);
    ylabelmine(temp2Lab,tagFont);
    
        set(gca,'Fontsize',axFont,'YColor','k') 
    
    fname=sprintf('ssfrTemp2_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'subdir',subDir,'fig','v'); end
    
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
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.','linewidth',lw,'color',ccc)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    set(gca,'YColor',ccc2)
    
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian(mmG),log10(yGal.yMedian(mmG)),'-','color',cc(1,:),'DisplayName','Galactic Gas','linewidth',lw);
    hold on
    h(2)=plot(yCGM.xMedian(mmC),log10(yCGM.yMedian(mmC)),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',lw);
    h(3)=plot(yOut.xMedian(mmO),log10(yOut.yMedian(mmO)),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',lw);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
    
     nTag=sprintf('%s, %s, %i',simDisplayName,mtag,sum(galaxyMask));
    xfac=0.6; yfac=1.025;
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',tagFont,'fontweight','bold')
    % annotations 
     plot(log10(ssfrThresh).*[1 1],yl,'k:','linewidth',lw);
    dx=diff(xl)/30;
    dy=diff(yl)/30;
        
    sfTagPos=[log10(ssfrThresh)+dx yl(2)-dy];
    
    text(sfTagPos(1),sfTagPos(2),sfTag,'color',[0 0 1],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    qTagPos =[log10(ssfrThresh)-6.6*dx yl(2)-dy];
    text(qTagPos(1),qTagPos(2),qTag,'color',[1 0 0],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    
    arrowLen=4*dx;
    sfAstart=[sfTagPos(1) sfTagPos(2)-dy];
    sfAend=[sfTagPos(1)+arrowLen sfTagPos(2)-dy];
    arrow(sfAstart,sfAend,'width',lw,'FaceColor',[0 0 1],'EdgeColor',[ 0 0 1])
    
    qAstart=[log10(ssfrThresh)-dx qTagPos(2)-dy];
    qAend=[qAstart(1)-arrowLen qTagPos(2)-dy];
    arrow(qAstart,qAend,'width',lw,'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
    
    grid
    
    %hl=legend(h);
    %set(hl,'Interpreter','latex','fontsize',tagFont,'location','SouthWest');
    
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,tagFont);
    ylabelmine(densLab,tagFont);
    
        set(gca,'Fontsize',axFont,'YColor','k') 
    fname=sprintf('ssfrDens_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'subdir',subDir,'fig','v'); end
    
     
    %% gas mass fraction density
    yl=[-3 0];
    yGal=mk_meanMedian_bin(xdata,fg(galaxyMask),'bins',bins);
    
    
    mmG=yGal.binCount>0;
    galPolY=log10([yGal.yQuarts(2,mmG) fliplr(yGal.yQuarts(3,mmG))]);
    galPolX=[yGal.xMedian(mmG) fliplr(yGal.xMedian(mmG))];
    galPolY(isinf(galPolY))=-10;
    
    
    figure('position',[1432 421 1000 750],'Color','w')
    yyaxis right
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.','linewidth',lw,'color',ccc)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    set(gca,'YColor',ccc2)
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian(mmG),log10(yGal.yMedian(mmG)),'-','color',cc(1,:),'DisplayName','Galactic Gas','linewidth',lw);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
     nTag=sprintf('%s, %s, %i',simDisplayName,mtag,sum(galaxyMask));
    xfac=0.6; yfac=1.025;
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',tagFont,'fontweight','bold')
    
    % annotations 
     plot(log10(ssfrThresh).*[1 1],yl,'k:','linewidth',lw);
    dx=diff(xl)/30;
    dy=diff(yl)/30;
        
    sfTagPos=[log10(ssfrThresh)+dx yl(2)-dy];
    
    text(sfTagPos(1),sfTagPos(2),sfTag,'color',[0 0 1],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    qTagPos =[log10(ssfrThresh)-6.6*dx yl(2)-dy];
    text(qTagPos(1),qTagPos(2),qTag,'color',[1 0 0],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    
    arrowLen=4*dx;
    sfAstart=[sfTagPos(1) sfTagPos(2)-dy];
    sfAend=[sfTagPos(1)+arrowLen sfTagPos(2)-dy];
    arrow(sfAstart,sfAend,'width',lw,'FaceColor',[0 0 1],'EdgeColor',[ 0 0 1])
    
    qAstart=[log10(ssfrThresh)-dx qTagPos(2)-dy];
    qAend=[qAstart(1)-arrowLen qTagPos(2)-dy];
    arrow(qAstart,qAend,'width',lw,'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
    
    grid
    
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,tagFont);
    ylabelmine(fgLab,tagFont);
    
        set(gca,'Fontsize',axFont,'YColor','k') 
    fname=sprintf('ssfrgasFrac_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'subdir',subDir,'fig','v'); end
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
    stairs(yGal.bins,yGal.binCount./sum(yGal.binCount),'-.','linewidth',lw,'color',ccc)
    ylabelmine('$N/N_\mathrm{Tot}$');
    ylim([0 1.05.*max(yGal.binCount./sum(yGal.binCount))])
    
    set(gca,'YColor',ccc2)
    yyaxis left
    h=[];
    h(1)=plot(yGal.xMedian(mmG),log10(yGal.yMedian(mmG)),'-','color',cc(1,:),'DisplayName','Galactic Gas','linewidth',lw);
    hold on
    h(2)=plot(yCGM.xMedian(mmC),log10(yCGM.yMedian(mmC)),'-','color',cc(2,:),'DisplayName','Inner CGM','linewidth',lw);
    h(3)=plot(yOut.xMedian(mmO),log10(yOut.yMedian(mmO)),'-','color',cc(3,:),'DisplayName','Outer CGM','linewidth',lw);
    
    patch(galPolX,galPolY,cc(1,:),'facealpha',0.4,'edgecolor','none')
    patch(cgmPolX,cgmPolY,cc(2,:),'facealpha',0.4,'edgecolor','none')
    patch(outPolX,outPolY,cc(3,:),'facealpha',0.4,'edgecolor','none')
     nTag=sprintf('%s, %s, %i',simDisplayName,mtag,sum(galaxyMask));
    xfac=0.6; yfac=1.025;
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),nTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',tagFont,'fontweight','bold')
    % annotations 
     plot(log10(ssfrThresh).*[1 1],yl,'k:','linewidth',lw);
    dx=diff(xl)/30;
    dy=diff(yl)/30;
        
    sfTagPos=[log10(ssfrThresh)+dx yl(2)-dy];
    
    text(sfTagPos(1),sfTagPos(2),sfTag,'color',[0 0 1],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    qTagPos =[log10(ssfrThresh)-6.6*dx yl(2)-dy];
    text(qTagPos(1),qTagPos(2),qTag,'color',[1 0 0],...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold');
    
    arrowLen=4*dx;
    sfAstart=[sfTagPos(1) sfTagPos(2)-dy];
    sfAend=[sfTagPos(1)+arrowLen sfTagPos(2)-dy];
    arrow(sfAstart,sfAend,'width',lw,'FaceColor',[0 0 1],'EdgeColor',[ 0 0 1])
    
    qAstart=[log10(ssfrThresh)-dx qTagPos(2)-dy];
    qAend=[qAstart(1)-arrowLen qTagPos(2)-dy];
    arrow(qAstart,qAend,'width',lw,'FaceColor',[1 0 0],'EdgeColor',[1 0 0])
    
    grid
    hl=legend(h);set(hl,'Interpreter','latex','fontsize',tagFont,'location','SouthWest');
    xlim(xl);ylim(yl);
    
    xlabelmine(ssfrLab,tagFont);
    ylabelmine(entLab,tagFont);
        set(gca,'Fontsize',axFont,'YColor','k') 
    
    fname=sprintf('ssfrEnt_%s',nameTag);
    if  printFlag; printout_fig(gcf,fname,'subdir',subDir,'fig','v'); end
    
    close all
end
