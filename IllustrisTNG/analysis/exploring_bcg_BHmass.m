
bp=illustris.set_env('100','nomount');



load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/BH_energyInjection_z0_TNG100.mat')
global simDisplayName
global DEFAULT_PRINTOUT_DIR

fofs=illustris.groupcat.loadHalos(bp,snap);
subs=illustris.groupcat.loadSubhalos(bp,snap);
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);

centralMask= subsInfo.isCentral(bhStruct.galMask);
m200=fofs.Group_M_Crit200(subsInfo.hostFof(bhStruct.galMask)+1).*illUnits.massUnit;

clustMask=m200>1e14;
mask=(clustMask & centralMask);

%%
centralBH=bhStruct.inGal.bhMassMax(clustMask & centralMask);


hf(1)=figure;

loglog(m200(mask),centralBH,'o','markersize',10)
grid
set(gca,'Fontsize',14)
%ylim([1e9,1.5e10])

xlabelmine('$M_{\mathrm{200,c}}\,[\mathrm{M_\mathrm{\odot}}]$',14)
ylabelmine('$M_{\mathrm{BH}}\,[\mathrm{M_\mathrm{\odot}}]$',14)
titlemine('BH Mass of TNG100 central galaxies in $M_\mathrm{Halo}>10^{14}\,M_{\mathrm{\odot}}$')


%name=sprintf('%s/figFiles/central_BHmass_clusterMass_%s.fig',DEFUALT_PRINTOUT_DIR,simDisplayName);
%savefig(hf,name,'compact')

%printout_fig(gcf,'central_BHmass_clusterMass_TNG100')
%

% hold on
% centralBH2=bhStruct.inGal.bhMassSum(clustMask & centralMask)
% centralBH=bhStruct.inGal.bhMassMax(clustMask & centralMask);
% loglog(m200(mask),centralBH2,'or','markersize',10)
% figure
% loglog(centralBH,centralBH2)
% loglog(centralBH,centralBH2,'.')
% hold on
% loglog([1e9 2e10],[1e9 2e10],'--k')
% loglog(centralBH,centralBH2,'o')
% centralBH./centralBH2-1


hf(2)=figure;

mgal=bhStruct.galMass(bhStruct.galMask);
loglog(mgal(mask),bhStruct.inGal.bhMassMax(mask),'o','markersize',10)
grid
set(gca,'Fontsize',14)
xlabelmine('$M_{\mathrm{*}}\,[\mathrm{M_\mathrm{\odot}}]$',14)
ylabelmine('$M_{\mathrm{BH}}\,[\mathrm{M_\mathrm{\odot}}]$',14)
titlemine('BH Mass of TNG100 central galaxies in $M_\mathrm{Halo}>10^{14}\,M_{\mathrm{\odot}}$')

 
    name=sprintf('%s/figFiles/central_BHmass_%s.fig',DEFAULT_PRINTOUT_DIR,simDisplayName);
    
    savefig(hf,name,'compact')

