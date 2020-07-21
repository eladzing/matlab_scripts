
bp=illustris.set_env(simname);


global simDisplayName

fname=['gasFraction_hostR200%s_' simDisplayName];
global illUnits
%global cosmoStruct

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

cc=brewermap(8,'Set1');


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
    text(xt,yt,str,'fontsize',20,'interpreter','latex')
end
grid
box on

% xfac=0.75; yfac=0.25;
% text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),tag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
% 'Interpreter','latex','fontsize',28,'fontweight','bold')
% 
zTag=sprintf('z=%3.2f',illUnits.zred);
xfac=0.05; yfac=0.075;
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


printout_fig(gcf,fname,'v')
