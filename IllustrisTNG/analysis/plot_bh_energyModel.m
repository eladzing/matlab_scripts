

global simDisplayName

if readFlag
    
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/cooling_times_z0_' simDisplayName '.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_z0_' simDisplayName '.mat'])
     
    fofs=illustris.groupcat.loadHalos(bp,99);
    subs=illustris.groupcat.loadSubhalos(bp,99);
end
units; % load general unit structure in cgs.
global illUnits

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.


%% define what we are plotting
% identify centrals
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

galMass=tCoolStruct.galMass(tCoolStruct.galMask);  % galaxy stellar mass


%% black hole stuff
bhQM=bhStruct.inGal.cumEngQM(tCoolStruct.galMask);
bhRM=bhStruct.inGal.cumEngRM(tCoolStruct.galMask);
bhMass=bhStruct.inGal.bhMassMax(tCoolStruct.galMask);

galaxyMask=centralMask & bhMass>0;
%% 
glm=galMass(galaxyMask);
bhm=bhMass(galaxyMask);
qm0=bhQM(galaxyMask).*illUnits.BHEnergyFactor; % in 1e53 ergs
km0=bhRM(galaxyMask).*illUnits.BHEnergyFactor; % in 1e53 ergs

qm=qm0;
km=km0;


qm(qm0==0)=1e-1;
km(km0==0)=1e-2.*10.^(0.5.*rand(1,sum(km0==0)));

xdataG=log10(glm);
xdataB=log10(bhm);
ydataQ=log10(qm);
ydataK=log10(km);
ydataR=log10(km./qm);

% do contours 
filt=fspecial('disk',4);
glmQcont=plot_population_contour(xdataG,ydataQ,'smooth',filt,'noplot');
glmKcont=plot_population_contour(xdataG,ydataK,'smooth',filt,'noplot');

bhmQcont=plot_population_contour(xdataB,ydataQ,'smooth',filt,'noplot');
bhmKcont=plot_population_contour(xdataB,ydataK,'smooth',filt,'noplot');

gmbmCont=plot_population_contour(xdataG,xdataB,'smooth',filt,'noplot');

bhRatCont=plot_population_contour(xdataB,ydataR,'smooth',filt,'noplot');

cmap=brewermap(4,'Set1');

lvl=8:15:100;

gmLab='$\log$ Stellar Mass $[\mathrm{M_\odot}]$';
bhLab='$\log$ BH Mass $[\mathrm{M_\odot}]$';
eLab='$\log$ Cum.\@ BH Output Energy  $[10^{53}\mathrm{erg}]$';
ratLab='$\log E_\mathrm{LAM}/E_\mathrm{HAM}$';
%kmLab='$\log$ Cum.\@ Energy  $[10^{53}\mathrm{erg}]$';




%%


%% 
figure
set(gcf,'position',[1432 421 1000 750],'Color','w')

h(1)=plot(xdataB,ydataQ,'o','color',cmap(2,:),'markerfacecolor',cmap(2,:),...
    'markersize',5,'DisplayName','HAM');
hold on 
h(2)=plot(xdataB,ydataK,'o','color',cmap(1,:),'markerfacecolor',cmap(1,:),...
    'markersize',5,'DisplayName','LAM');

contour(bhmQcont.xx,bhmQcont.yy,bhmQcont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

contour(bhmKcont.xx,bhmKcont.yy,bhmKcont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

grid

hl=legend(h);
set(hl,'Interpreter','Latex','Fontsize',14,'location','SouthEast')

set(gca,'Fontsize',14)
xl=[6 10.5];
yl=[-2.1 10.2];

xlim(xl)
ylim(yl)

xlabelmine(bhLab);
ylabelmine(eLab);

fname=sprintf('bhModel_bhMass_%s',simDisplayName);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end
%%

hf=figure;
set(gcf,'position',[1432 421 1000 750],'Color','w')
axes1 = axes('Parent',hf);
hold(axes1,'on');

plot(xdataG,ydataQ,'o','color',cmap(2,:),'markerfacecolor',cmap(2,:),...
    'markersize',5,'DisplayName','HAM');
plot(xdataG,ydataK,'o','color',cmap(1,:),'markerfacecolor',cmap(1,:),...
    'markersize',5,'DisplayName','LAM');

contour(glmQcont.xx,glmQcont.yy,glmQcont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

contour(glmKcont.xx,glmKcont.yy,glmKcont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

grid

% hl=legend(h);
% set(hl,'Interpreter','Latex','Fontsize',14,'location','SouthEast')

xl=[9 12.5];
yl=[-2.1 10.2];
set(gca,'Fontsize',14)
xlim(xl)
ylim(yl)

xlabelmine(gmLab,16);
ylabelmine(eLab,16);

box(axes1,'on');
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14);

axes2 = axes('Parent',hf,...
    'Position',[0.568799298860649 0.194472049689441 0.313413672217354 0.380683229813665]);
hold(axes2,'on');


plot(xdataG,xdataB,'o','color',cmap(2,:),'markerfacecolor',cmap(2,:),...
    'markersize',5,'DisplayName','HAM');
hold on 

contour(gmbmCont.xx,gmbmCont.yy,gmbmCont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

box(axes2,'on');
grid(axes2,'on');
% Set the remaining axes properties
set(axes2,'FontSize',14);


yl=[6 10.5];
xl=[9 12.5];

xlim(xl)
ylim(yl)


xlabelmine(gmLab);
ylabelmine(bhLab);

fname=sprintf('bhModel_stellarMass_%s',simDisplayName);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end


%%

figure
set(gcf,'position',[1432 421 1000 750],'Color','w')

h(1)=plot(xdataB,ydataR,'o','color',cmap(2,:),'markerfacecolor',cmap(2,:),...
    'markersize',5,'DisplayName','HAM');
hold on 
% h(2)=plot(xdataB,ydataK,'o','color',cmap(1,:),'markerfacecolor',cmap(1,:),...
%     'markersize',5,'DisplayName','LAM');
% 
 contour(bhRatCont.xx,bhRatCont.yy,bhRatCont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

% contour(bhmKcont.xx,bhmKcont.yy,bhmKcont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
%     'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

grid

hl=legend(h);
set(hl,'Interpreter','Latex','Fontsize',14,'location','SouthEast')

set(gca,'Fontsize',14)
xl=[6 10.5];
yl=[-2.1 10.2];

xlim(xl)
ylim(yl)

xlabelmine(bhLab);
ylabelmine(ratLab);

fname=sprintf('bhModel_bhMass_%s',simDisplayName);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end




    
    
    
    
% 
% loglog(glm,qm,'.',glm,km,'.')
% ylim([5e-3 1e10])
% 
% 
% 
% figure(2)
% loglog(bhm,qm,'.',bhm,km,'.')
% ylim([5e-3 1e10])
% 
% figure(3)
% loglog(glm,bhm,'.')
% 
% 
% 
% 
% 




% %% plot number density plot 
% 
% 
% xyK=cat(1,xdata,ydataK);
% xyQ=cat(1,xdata,ydataQ);
% 
% %cmapN=brewermap(256,'*YlOrRd');
% cmapN=brewermap(256,'*Spectral');
% 
% [~,dd]=knnsearch(xyK',xyK','k',3);
% ddK=1./dd(:,end).^2;
% [~,dd]=knnsearch(xyQ',xyQ','k',3);
% ddQ=1./dd(:,end).^2;
% figure
% scatter(xdata,ydataQ,2,log10(ddQ))
% hold on
% scatter(xdata,ydataK,2,log10(ddK))
% caxis([1 5])
% colormap(cmapN)
% colorbar
% 
% 
% 
% %% plot contour plots 
% 
% 
% cmapQ=brewermap(256,'*YlGnBu');
% 
% % stellar mass 
% c1=glmQcont;
% msk=c1.popContour==100;
% c1.popContour(msk)=NaN;
% 
% c2=glmKcont;
% msk=c2.popContour==100;
% c2.popContour(msk)=NaN;
% 
% 
% xl=[9 12.5];
% yl=[-2 10];
% 
% figure 
% set(gcf,'position',[1432 421 1000 750],'Color','w')
% 
% ax0=axes; 
% xlim(xl);ylim(yl);
% 
% contourf(c1.xx,c1.yy,c1.popContour,'ShowText','off','LineColor',cmap(2,:),...
%     'LevelList',0:20:90,'Fill','on','linestyle','-','linewidth',1.5);
% colormap(cmapQ);
% 
% ax2=axes; 
% contourf(c2.xx,c2.yy,c2.popContour,'ShowText','off','LineColor',cmap(1,:),...
%     'LevelList',[0:20:90],'Fill','on','linestyle','-','linewidth',1.5);
% colormap(cmapK);
% 
% patches = findobj(ch,'-property','AlphaData');
% for ph = patches
%  set(ph,'AlphaData', 0.1 * get(ph,'AlphaData'));
% end
% 
% xlim(xl);ylim(yl);
% 
% set(ax2,'position',get(ax0,'position'));
% linkaxes([ax0,ax2])
% ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
% ax1.Visible = 'off';ax1.XTick = [];ax1.YTick = [];
% 
% 
% 
% 
% xlabelmine(bhLab,20);
% ylabelmine(entLab,20);
% set(gca,'fontsize',18)
% 
