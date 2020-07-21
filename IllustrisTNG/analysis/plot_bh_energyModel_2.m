
global DRACOFLAG
global simDisplayName

if readFlag
    
    global DEFAULT_MATFILE_DIR
%     load([DEFAULT_MATFILE_DIR '/cooling_times_z0_' simDisplayName '.mat'])
%     load([DEFAULT_MATFILE_DIR '/BH_energyInjection_z0_' simDisplayName '.mat'])
%      
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_snp' num2str(snap) '_' simDisplayName '.mat'])
     load([DEFAULT_MATFILE_DIR '/gasProperties_snp' num2str(snap) '_' simDisplayName '.mat'])
 
     if DRACOFLAG
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
    else
        loadFofSubTNG100
    end
    
end
units; % load general unit structure in cgs.
global illUnits

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.

massThresh=10^9; % threshold for *stellar* mass


%% define what we are plotting
% identify centrals
%centralMask= subsInfo.isCentral(tCoolStruct.galMask);

galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'centrals');
spbMask = mk_splashback_mask('time',5,'both',0.1);

galMass=tCoolStruct.galMass;  % galaxy stellar mass
bhMass=bhStruct.inGal.bhMassMax;

galaxyMask=galMask & ~spbMask & bhMass>0;

%% black hole stuff
bhQM=bhStruct.inGal.cumEngQM(galaxyMask).*illUnits.BHEnergyFactor; % in 1e53 ergs;
bhRM=bhStruct.inGal.cumEngRM(galaxyMask).*illUnits.BHEnergyFactor; % in 1e53 ergs;
%bhMass=bhStruct.inGal.bhMassMax(galaxyMask);

%% 
glm=galMass(galaxyMask);
bhm=bhMass(galaxyMask);
qm0=bhQM;  %(galaxyMask).*illUnits.BHEnergyFactor; % in 1e53 ergs
km0=bhRM;  %(galaxyMask).*illUnits.BHEnergyFactor; % in 1e53 ergs

qm=qm0;
km=km0;


qm(qm0==0)=1e-1;
km(km0==0)=1e-2.*10.^(0.5.*rand(1,sum(km0==0)));
bhRat=km0./qm0;
bhRat(km0==0 | qm0==0)=1e-8.*10.^(0.5.*rand(1,sum(km0==0 | qm0==0  )));

xdataG=log10(glm);
xdataB=log10(bhm);
ydataQ=log10(qm);
ydataK=log10(km);
ydataR=log10(bhRat);


% do contours 
filt=fspecial('disk',4);
glmQcont=plot_population_contour(xdataG,ydataQ,'smooth',filt,'noplot');
glmKcont=plot_population_contour(xdataG,ydataK,'smooth',filt,'noplot');

bhmQcont=plot_population_contour(xdataB,ydataQ,'smooth',filt,'noplot');
bhmKcont=plot_population_contour(xdataB,ydataK,'smooth',filt,'noplot');

gmbmCont=plot_population_contour(xdataG,xdataB,'smooth',filt,'noplot');

bhRatCont=plot_population_contour(xdataB,ydataR,'smooth',filt,'noplot');

%cmap=brewermap(4,'Set1');

lvl=8:15:100;

gmLab='$\log$ Stellar Mass $[\mathrm{M_\odot}]$';
bhLab='$\log$ BH Mass $[\mathrm{M_\odot}]$';
eLab='$\log$ Cum.\@ BH Output Energy  $[10^{53}\mathrm{erg}]$';
ratLab='$\log E_\mathrm{Kinetic}/E_\mathrm{Thermal}$';
%kmLab='$\log$ Cum.\@ Energy  $[10^{53}\mathrm{erg}]$';




%% population color maps 
xlG=[9 12.5];
xlB=[6 10.5];
yl=[-2.1 10.2];
ylR=[-8 0];

bs=100;

[dataGQ, ~, ~,~]= histogram2d(xdataG,ydataQ,ones(size(xdataG)),...
    'xxlim',xlG,'yylim',yl,'bins',bs);
[dataGK, ~, ~,~]= histogram2d(xdataG,ydataK,ones(size(xdataG)),...
    'xxlim',xlG,'yylim',yl,'bins',bs);

[dataBQ, ~, ~,~]= histogram2d(xdataB,ydataQ,ones(size(xdataB)),...
    'xxlim',xlB,'yylim',yl,'bins',bs);
[dataBK, ~, ~,~]= histogram2d(xdataB,ydataK,ones(size(xdataB)),...
    'xxlim',xlB,'yylim',yl,'bins',bs);

[dataGB, ~, ~,~]= histogram2d(xdataG,xdataB,ones(size(xdataG)),...
    'xxlim',xlG,'yylim',xlB,'bins',bs);
[dataBR, ~, ~,~]= histogram2d(xdataB,ydataR,ones(size(xdataB)),...
    'xxlim',xlB,'yylim',ylR,'bins',bs);

dataGQ=(dataGQ(:,:,1));
dataGK=(dataGK(:,:,1));

dataBQ=(dataBQ(:,:,1));
dataBK=(dataBK(:,:,1));

dataGB=(dataGB(:,:,1));
dataBR=(dataBR(:,:,1));

dGQ=min(dataGQ,1);
dGK=min(dataGK,1);
dGK(dGK==1)=-1;


dBQ=min(dataBQ,1);
dBK=min(dataBK,1);
dBK(dBK==1)=-1;

dGB=min(dataGB,1);dGB(dGB==1)=1.66; %1.2
dBR=min(dataBR,1);dBR(dBR==1)=1.66; %1.2





concol=[0,0,0]; %brewermap(1,'Dark2');
%% plotting 

if ~DRACOFLAG
    tagFont=18;
    axFont=18;
    bigTagFont=20;
    lw=3;
else
    tagFont=30;
    axFont=30;
    bigTagFont=34;
    lw=5;
end


cmap=brewermap(12,'RdBu');
cmap(1,:)=[1 1 1];

hf=figure;
set(gcf,'position',[1432 421 1000 750],'Color','w')
axes1 = axes('Parent',hf);
%hold(axes1,'on');
h=[];
data=dGQ+dGK;
data(data==0)=-1;
data(dGQ==0 & dGK==0)=-10;
imagesc(xlG,yl,data)
caxis([-1.8 1.8])
colormap(cmap)
hold on
contour(glmQcont.xx,glmQcont.yy,glmQcont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',1);

contour(glmKcont.xx,glmKcont.yy,glmKcont.popContour,'ShowText','off','LineColor',concol,...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',1);

h(1)=plot(xlG(1)-1,yl(1)-1,'s','color',cmap(10,:),'markerfacecolor',cmap(10,:),...
'markersize',15,'DisplayName','Thermal Mode');
h(2)=plot(xlG(1)-1,yl(1)-1,'s','color',cmap(3,:),'markerfacecolor',cmap(3,:),...
'markersize',15,'DisplayName','Kinetic Mode');

hl=legend(h);
set(hl,'Interpreter','Latex','Fontsize',tagFont,'location','NorthWest')

set(gca,'Ydir','normal','Fontsize',axFont)
xlim(xlG)
ylim(yl)
grid

xlabelmine(gmLab,tagFont);
ylabelmine(eLab,tagFont);

box(axes1,'on');


%grid(axes1,'on');

% Set the remaining axes properties
set(axes1,'FontSize',axFont);

axes2 = axes('Parent',hf,...
    'Position',[0.568799298860649 0.194472049689441 0.313413672217354 0.380683229813665]);
%hold(axes2,'on');


imagesc(xlG,xlB,dGB)
set(axes2,'Ydir','normal','FontSize',tagFont);
hold on 

contour(gmbmCont.xx,gmbmCont.yy,gmbmCont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);
caxis([0.6 2.4])
box(axes2,'on');
grid(axes2,'on');
% Set the remaining axes properties



xlim(xlG)
ylim(xlB)

%xlabelmine(gmLab);
hlab=ylabelmine(bhLab,tagFont);
set(hlab,'backgroundcolor',[1 1 1]);


fname=sprintf('bhModel_stellarMass_%s',simDisplayName);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end


%% 

hf=figure;
set(gcf,'position',[1432 421 1000 750],'Color','w')
axes1 = axes('Parent',hf);
%hold(axes1,'on');
h=[];
data=dBQ+dBK;
data(data==0)=-1;
data(dBQ==0 & dBK==0)=-10;

imagesc(xlB,yl,data)
caxis([-1.8 1.8])
colormap(cmap)
hold on
contour(bhmQcont.xx,bhmQcont.yy,bhmQcont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

contour(bhmKcont.xx,bhmKcont.yy,bhmKcont.popContour,'ShowText','off','LineColor',concol,...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);

h(1)=plot(xlB(1)-1,yl(1)-1,'s','color',cmap(10,:),'markerfacecolor',cmap(10,:),...
'markersize',15,'DisplayName','Thermal Mode');
h(2)=plot(xlB(1)-1,yl(1)-1,'s','color',cmap(3,:),'markerfacecolor',cmap(3,:),...
'markersize',15,'DisplayName','Kinetic Mode');

hl=legend(h);
set(hl,'Interpreter','Latex','Fontsize',tagFont,'location','NorthWest')

set(gca,'Ydir','normal','Fontsize',axFont)
xlim(xlB)
ylim(yl)
grid

xlabelmine(bhLab,tagFont);
ylabelmine(eLab,tagFont);

box(axes1,'on');


%grid(axes1,'on');

% Set the remaining axes properties
set(axes1,'FontSize',axFont);

axes2 = axes('Parent',hf,...
    'Position',[0.568799298860649 0.194472049689441 0.313413672217354 0.380683229813665]);
%hold(axes2,'on');


imagesc(xlB,ylR,dBR)
set(axes2,'Ydir','normal','FontSize',axFont);
hold on 

contour(bhRatCont.xx,bhRatCont.yy,bhRatCont.popContour,'ShowText','off','LineColor',[0 0 0 ],...
    'LevelList',lvl,'Fill','off','linestyle','-','linewidth',0.5);
caxis([0.6 2.4])
box(axes2,'on');
grid(axes2,'on');
% Set the remaining axes properties



xlim(xlB)
ylim(ylR)

hlab=ylabelmine(ratLab,tagFont);
set(hlab,'backgroundcolor',[1 1 1]);
%xlabelmine(bhLab);

fname=sprintf('bhModel_bhMass_%s',simDisplayName);
    if  printFlag; printout_fig(gcf,fname,'fig','v'); end












