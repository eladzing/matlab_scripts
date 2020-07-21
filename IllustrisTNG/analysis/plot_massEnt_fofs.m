
bp=illustris.set_env(simname);

global simDisplayName


global illUnits
global cosmoStruct

if readFlag
    
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/gasProperties_Fofs_snp' num2str(snap) '_' simDisplayName '.mat'])

   fofs=illustris.groupcat.loadHalos(bp,99);
   subs=illustris.groupcat.loadSubhalos(bp,99);

end


subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
hostMass=fofs.Group_M_Crit200(subsInfo.hostFof(gasPropsFofs.galMask)+1).*illUnits.massUnit;
galMass=gasPropsFofs.galMass(gasPropsFofs.galMask);
gmass=gasPropsFofs.inR200.gasMass(gasPropsFofs.galMask)+gasPropsFofs.inR200.sfrMass(gasPropsFofs.galMask);
ent=gasPropsFofs.inR200.meanEntMW(gasPropsFofs.galMask);

hostm=mk_meanMedian_bin(log10(galMass),log10(hostMass),'bins',9:0.5:12.5);
enth1=mk_meanMedian_bin(log10(hostMass),log10(ent),'bins',11:0.25:15);
entg1=mk_meanMedian_bin(log10(galMass),log10(ent),'bins',9:0.5:12.5);


% mh=[min(hostMass) max(hostMass)];
mh=[3e10 1e15];
kk=mass_entropy_relation(mh,'cosmo',cosmoStruct);
[gmin,im]=min(galMass);
[gmax,ix]=max(galMass);
kkg=mass_entropy_relation(hostMass([im ix]),'cosmo',cosmoStruct);
kkg2=mass_entropy_relation(hostMass,'cosmo',cosmoStruct);

entg2=mk_meanMedian_bin(log10(galMass),log10(kkg2),'bins',9:0.5:12.5);


%% fiducial
snap35=4;

bp1=illustris.set_env(35,'0000');

load([DEFAULT_MATFILE_DIR '/methods/gasProperties_fofs_snp4_' simDisplayName '.mat'])

fofs35=illustris.groupcat.loadHalos(bp1,snap35);
subs35=illustris.groupcat.loadSubhalos(bp1,snap35);
subsInfo35 = illustris.infrastructure.build_sub_fof_connection(subs35,fofs35);

hostMassMeth.fid=fofs35.Group_M_Crit200(subsInfo35.hostFof(gasPropsFofs.galMask)+1).*illUnits.massUnit;
galMassMeth.fid=gasPropsFofs.galMass(gasPropsFofs.galMask);
gmassMeth.fid=gasPropsFofs.inR200.gasMass(gasPropsFofs.galMask)+gasPropsFofs.inR200.sfrMass(gasPropsFofs.galMask);
entMeth.fid=gasPropsFofs.inR200.meanEntMW(gasPropsFofs.galMask);



% hostm=mk_meanMedian_bin(log10(galMass),log10(hostMass),'bins',9:0.5:12.5);
% ent1_35=mk_meanMedian_bin(log10(hostMass35),log10(ent35),'bins',11:0.25:15);
%     

%% no BH 
bp1=illustris.set_env(35,'2201');

load([DEFAULT_MATFILE_DIR '/methods/gasProperties_fofs_snp4_' simDisplayName '.mat'])

fofs35=illustris.groupcat.loadHalos(bp1,snap35);
subs35=illustris.groupcat.loadSubhalos(bp1,snap35);
subsInfo35 = illustris.infrastructure.build_sub_fof_connection(subs35,fofs35);

hostMassMeth.noBH=fofs35.Group_M_Crit200(subsInfo35.hostFof(gasPropsFofs.galMask)+1).*illUnits.massUnit;
galMassMeth.noBH=gasPropsFofs.galMass(gasPropsFofs.galMask);
gmassMeth.noBH=gasPropsFofs.inR200.gasMass(gasPropsFofs.galMask)+gasPropsFofs.inR200.sfrMass(gasPropsFofs.galMask);
entMeth.noBH=gasPropsFofs.inR200.meanEntMW(gasPropsFofs.galMask);
    
%% no Lam

bp1=illustris.set_env(35,'3000');

load([DEFAULT_MATFILE_DIR '/methods/gasProperties_fofs_snp4_' simDisplayName '.mat'])

fofs35=illustris.groupcat.loadHalos(bp1,snap35);
subs35=illustris.groupcat.loadSubhalos(bp1,snap35);
subsInfo35 = illustris.infrastructure.build_sub_fof_connection(subs35,fofs35);

hostMassMeth.noLam=fofs35.Group_M_Crit200(subsInfo35.hostFof(gasPropsFofs.galMask)+1).*illUnits.massUnit;
galMassMeth.noLam=gasPropsFofs.galMass(gasPropsFofs.galMask);
gmassMeth.noLam=gasPropsFofs.inR200.gasMass(gasPropsFofs.galMask)+gasPropsFofs.inR200.sfrMass(gasPropsFofs.galMask);
entMeth.noLam=gasPropsFofs.inR200.meanEntMW(gasPropsFofs.galMask);
    

%% adiabatic
bp1=illustris.set_env(35,'0030');

load([DEFAULT_MATFILE_DIR '/methods/gasProperties_fofs_snp4_' simDisplayName '.mat'])

fofs35=illustris.groupcat.loadHalos(bp1,snap35);
subs35=illustris.groupcat.loadSubhalos(bp1,snap35);
subsInfo35 = illustris.infrastructure.build_sub_fof_connection(subs35,fofs35);

hostMassMeth.adiabat=fofs35.Group_M_Crit200(subsInfo35.hostFof(gasPropsFofs.galMask)+1).*illUnits.massUnit;
galMassMeth.adiabat=gasPropsFofs.galMass(gasPropsFofs.galMask);
gmassMeth.adiabat=gasPropsFofs.inR200.gasMass(gasPropsFofs.galMask)+gasPropsFofs.inR200.sfrMass(gasPropsFofs.galMask);
entMeth.adiabat=gasPropsFofs.inR200.meanEntMW(gasPropsFofs.galMask);
    

%% plot 

entgFid=mk_meanMedian_bin(log10(galMassMeth.fid),log10(entMeth.fid),'bins',9:0.5:12.5);
entgNobh=mk_meanMedian_bin(log10(galMassMeth.noBH),log10(entMeth.noBH),'bins',9:0.5:12.5);
entgNolam=mk_meanMedian_bin(log10(galMassMeth.noLam),log10(entMeth.noLam),'bins',9:0.5:12.5);
entgAdiabat=mk_meanMedian_bin(log10(galMassMeth.adiabat),log10(entMeth.adiabat),'bins',9:0.5:12.5);

enthFid=mk_meanMedian_bin(log10(hostMassMeth.fid),log10(entMeth.fid),'bins',11:0.25:15);
enthNobh=mk_meanMedian_bin(log10(hostMassMeth.noBH),log10(entMeth.noBH),'bins',11:0.25:15);
enthNolam=mk_meanMedian_bin(log10(hostMassMeth.noLam),log10(entMeth.noLam),'bins',11:0.25:15);
enthAdiabat=mk_meanMedian_bin(log10(hostMassMeth.adiabat),log10(entMeth.adiabat),'bins',9:0.5:12.5);

cc=brewermap(8,'Set1');
%     host

% 
% 
% plot(log10(hostMass),log10(ent),'.');
% hold on
% plot(ent1.xMedian,ent1.yMedian)
% 
% plot(log10(mh),log10(kk))


figure('color','w')
h(1)=errorbar(entg1.xMedian,entg1.yMedian,...
    abs(entg1.yQuarts(1,:)-entg1.yMedian),...
    abs(entg1.yQuarts(4,:)-entg1.yMedian),...
    '-','DisplayName','TNG100','linewidth',1.5,'color',cc(1,:));
hold on
h(2)=errorbar(entgNobh.xMedian,entgNobh.yMedian,...
    abs(entgNobh.yQuarts(1,:)-entgNobh.yMedian),...
    abs(entgNobh.yQuarts(4,:)-entgNobh.yMedian),...
    '-r','DisplayName','no BH','linewidth',1.5,'color',cc(2,:));
h(3)=errorbar(entgNolam.xMedian,entgNolam.yMedian,...
    abs(entgNolam.yQuarts(1,:)-entgNolam.yMedian),...
    abs(entgNolam.yQuarts(4,:)-entgNolam.yMedian),...
    '-g','DisplayName','no LAM','linewidth',1.5,'color',cc(3,:));
h(4)=errorbar(entg2.xMedian,entg2.yMedian,...
    abs(entg2.yQuarts(1,:)-entg2.yMedian),...
    abs(entg2.yQuarts(4,:)-entg2.yMedian),...
    '-','DisplayName','Adiabatic','linewidth',1.5,'color','k');
%plot(log10(galMass),log10(kkg2),'.k','linewidth',2,...
%    'DisplayName','Adiabatic');


grid 
hl=legend(h);
set(hl,'fontsize',14,'interpreter','latex','location','SouthEast')

xlabelmine('$\log M_\mathrm{gal}\,[\mathrm{M_\odot}]$');
ylabelmine('$\log K\,[\mathrm{KeV\,cm^2}]$');

set(gca,'fontsize',14)
    
xlim([9 12.5])
printout_fig(gcf,'gMass_fofEnt_snp99_TNG100')

figure('color','w')

h(1)=errorbar(enth1.xMedian,enth1.yMedian,...
    abs(enth1.yQuarts(1,:)-enth1.yMedian),...
    abs(enth1.yQuarts(4,:)-enth1.yMedian),...
    '-b','DisplayName','TNG100','linewidth',1.5,'color',cc(1,:));
hold on
h(2)=errorbar(enthNobh.xMedian,enthNobh.yMedian,...
    abs(enthNobh.yQuarts(1,:)-enthNobh.yMedian),...
    abs(enthNobh.yQuarts(4,:)-enthNobh.yMedian),...
    '-r','DisplayName','no BH','linewidth',1.5,'color',cc(2,:));
h(3)=errorbar(enthNolam.xMedian,enthNolam.yMedian,...
    abs(enthNolam.yQuarts(1,:)-enthNolam.yMedian),...
    abs(enthNolam.yQuarts(4,:)-enthNolam.yMedian),...
    '-g','DisplayName','no LAM','linewidth',1.5,'color',cc(3,:));
h(4)=plot(log10(mh),log10(kk),'--k','linewidth',2,...
    'DisplayName','Adiabatic');


grid 
hl=legend(h);
set(hl,'fontsize',14,'interpreter','latex','location','SouthEast')

xlabelmine('$\log M_\mathrm{host}\,[\mathrm{M_\odot}]$');
ylabelmine('$\log K\,[\mathrm{KeV\,cm^2}]$');

set(gca,'fontsize',14)
    
xlim([11 15])
printout_fig(gcf,'hMass_fofEnt_snp99_TNG100')







% h(1)=errorbar(ent1.xMedian,ent1.yMedian,)
% hold on
% plot(entFid.xMedian,entFid.yMedian)
% plot(entNobh.xMedian,entNobh.yMedian)
% plot(entNolam.xMedian,entNolam.yMedian)
% 
% 
% plot(log10(mh),log10(kk))
% 
% 
% 
% plot(log10(hostMass35),log10(ent35),'.');
% hold on
% 
% 








filt=fspecial('disk',6);
xdata=log10(galMass);
ydata=log10(gmass./hostMass);
%obal popCont
popCont=plot_population_contour(xdata,ydata,'smooth',filt,'noplot');
yl=[-3 -0.5];
xl=[9 12.5];

minLev=7;
splitParam=30;
tre=tree2d(xdata,ydata,[xl(1) yl(1)],diff(xl) ,diff(yl), splitParam,minLev);

cax=[1 3]; %[min(cc) max(cc)];

%% plotting 
figure
set(gcf,'position',[1432 421 1000 750],'Color','w')


ax0=axes; 
xlim(xl);ylim(yl);
len=length(hostm.binEdges)-1;
cm=brewermap(len+1,'RdPu');

yp=[yl(1) yl(1) yl(2) yl(2)];
yt=yl(2)-0.05*diff(yl);

for i=1:len
    str=sprintf('$10^{%3.1f}$',hostm.yMedian(i));
    xt=hostm.binEdges(i)+0.3.*diff(hostm.binEdges(i:i+1));
    xp=cat(2,hostm.binEdges(i:i+1),fliplr(hostm.binEdges(i:i+1)));
    patch(xp,yp,cm(i,:),'facealpha',0.2,'edgecolor','k','edgealpha',0.4)
    hold on
    text(xt,yt,str,'fontsize',16,'interpreter','latex')
end
grid
box on

% xfac=0.75; yfac=0.25;
% text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),tag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
% 'Interpreter','latex','fontsize',28,'fontweight','bold')
% 
zTag=sprintf('z=%3.2f',illUnits.zred);
xfac=0.86; yfac=0.075;
text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),zTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
'Interpreter','latex','fontsize',24,'fontweight','bold')

mstarLab='log Stellar Mass $[\mathrm{M_\odot}]$';
fgLab='$\log M_\mathrm{gas}/M_\mathrm{200,c}$';

xlabelmine(mstarLab,20);
ylabelmine(fgLab,20);
set(gca,'fontsize',20)

ax1=axes;

cmap=brewermap(256,'YlOrRd');

celVal=points2tree(log10(ent),tre,'median');
plot2dTreeMap(celVal,tre,'fig',gcf,'cmap',cmap,'minmax',cax);

hold on

contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');

% qCol=[0.89,0.1,0.11];
% [~,h(1)]=contour(qCont.xx,qCont.yy,qCont.popContour,'ShowText','off',...
%     'LineColor',qCol,'linewidth',8,...
%     'LevelList',90,'Fill','off','linestyle','--',...
%    'DisplayName','Quenched');
% 
% sCol=[0.1,0.1,0.89];
% [~,h(2)]=contour(sfCont.xx,sfCont.yy,sfCont.popContour,'ShowText','off',...
%     'LineColor',sCol,'linewidth',8,...
%     'LevelList',90,'Fill','off','linestyle','--',...
%     'DisplayName','Star-forming');

% 
% xfac=0.1; yfac=0.85;
% text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),'Quenched',...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
% 'Interpreter','latex','fontsize',32,'fontweight','bold','color',qCol)
%  
% xfac=0.5; yfac=0.1;
% text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),'Star-forming',...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
% 'Interpreter','latex','fontsize',32,'fontweight','bold','color',sCol)
%  

barLab='$\log K \,[\mathrm{KeV\,cm^2}]$';

hb=colorbar;barTitle(hb,barLab,'fontsize',20)
set(hb,'fontsize',20);
colormap(cmap);caxis(cax);
xlim(xl);ylim(yl);

set(ax1,'position',get(ax0,'position'));
linkaxes([ax0,ax1])  %linkaxes([ax0,ax1,ax2])
%ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
ax1.Visible = 'off';ax1.XTick = [];ax1.YTick = [];

fname=['gasFraction_hostR200_' simDisplayName];
printout_fig(gcf,fname,'v')
