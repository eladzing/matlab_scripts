snap=99;
global simDisplayName
global illUnits
illustris.utils.set_illUnits(snap);

if readFlag
      
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
end





massThresh=10^9; % threshold for *stellar* mass

subsInfo=illustris.infrastructure.build_sub_fof_connection(subs,fofs);
spbMask = mk_splashback_mask('time',5,'both',0.1);

gMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');

galMask=gMask & ~spbMask;


sMassAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
gMmassAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf

galMass=sMassAllGals(galMask);
fgs=gMmassAllGals(galMask)./sMassAllGals(galMask);

mfg=fgs>0;
xdata=log10(galMass(mfg));
ydata=log10(fgs(mfg));
filt=fspecial('disk',6);
pCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
med = mk_meanMedian_bin(xdata,ydata,'nb',50);

figure
contour(pCont.xx,pCont.yy,pCont.popContour,'ShowText','off','LineColor',[0 0 0],...
                'LevelList',20:20:100,'Fill','off','linestyle','-');
  hold on
 plot(med.xMedian,med.yMedian,'-r','linewidth',1.6)
 plot(med.xMedian,med.yQuarts(2,:),'--r','linewidth',1.2)
 plot(med.xMedian,med.yQuarts(3,:),'--r','linewidth',1.2)
 
 grid
 set(gca,'fontsize',14);
 xlabelmine('log stellar mass');
 ylabelmine('log gas-to-stellar mass ratio');
 titlemine(simDisplayName);

printout_fig(gcf,['fgs_stellarMass_' simDisplayName],'v')

save