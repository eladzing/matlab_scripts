%% plot galaxy size vs stellar mass 

snap=99; %z=0
bp=illustris.set_env(simName);
global illUnits
if ~exist('readFlag','var')
    readFlag=true;
end

if readFlag
    fprintf(' *** Reading data *** \n');
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    readFlag=false;
    
    global simDisplayName
    global DEFAULT_MATFILE_DIR
    load([DEFAULT_MATFILE_DIR '/cooling_times_z0_' simDisplayName '.mat'])
    load([DEFAULT_MATFILE_DIR '/BH_energyInjection_z0_' simDisplayName '.mat'])
    
    
end


massThresh=10^9; % threshold for *stellar* mass


%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
units; % load general unit structure in cgs.

%% select galaxies
gasTemp=tCoolStruct.inGal.meanTempMW(1,:);


% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'centrals');




rhalfStar=subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,:); % stellar half mass radius
rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,:); % gas half mass radius


galMask=galMask & gasTemp>0 & rhalfGas>0 & rhalfStar>0;


stellarMass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galMask).*illUnits.massUnit; % stellar mass within 2*rhalf
hostMass=fofs.Group_M_Crit200(subsInfo.hostFof(galMask)+1).*illUnits.massUnit; % host mass 
hostRad=fofs.Group_R_Crit200(subsInfo.hostFof(galMask)+1).*illUnits.lengthUnit; % host mass 




rhalfStar=rhalfStar(galMask);
rhalfGas=rhalfGas(galMask);


rStar=2.*rhalfStar;


hostMed=mk_meanMedian_bin(log10(stellarMass),log10(hostMass),'bins',9:0.5:12.5);
hostMedR=mk_meanMedian_bin(log10(stellarMass),(hostRad),'bins',9:0.5:12.5);

%% stellar mass 


xdata=log10(stellarMass);

y1data=log10(rStar);
y2data=log10(rhalfGas);


filt=fspecial('disk',8);
lvl=8:15:100;
popCont1=plot_population_contour(xdata,y1data,'smooth',filt,'noplot');
popCont2=plot_population_contour(xdata,y2data,'smooth',filt,'noplot');

xl=[9 12.5];
yl=[0.4 3.2];


bs=100;
[dataS, ~, ~,~]= histogram2d(xdata,y1data,ones(size(xdata)),...
    'xxlim',xl,'yylim',yl,'bins',bs);
[dataG, ~, ~,~]= histogram2d(xdata,y2data,ones(size(xdata)),...
    'xxlim',xl,'yylim',yl,'bins',bs);

dataG=(dataG(:,:,1));
dataS=(dataS(:,:,1));
dG=min(dataG,1);
dS=min(dataS,1);
dS(dS==1)=-1;


mstarLab='$\log\mathrm{Stellar Mass}\,[\mathrm{M_\odot}]$';
ylab='$\log \mathrm{Radius}\,[\mathrm{kpc}]$';

concol=[0,0,0];
%% host mass 

yp=[yl(1) yl(1) yl(2) yl(2)];
yt=yl(2)-0.05*diff(yl);

len=length(hostMedr.binEdges)-1;
cm=brewermap(len+1,'RdPu');

cmap=brewermap(12,'RdBu');
cmap(1,:)=[1 1 1];

%% plot 

figure('position',[186   516   766   649],'Color','w');
ax0=axes; 


xlim(xl);ylim(yl);


for i=1:len
 str=sprintf('$%4.3f$',hostMedR.yMedian(i)/1e3);
    xt=hostMedR.binEdges(i)+0.25.*diff(hostMedR.binEdges(i:i+1));
    xp=cat(2,hostMedR.binEdges(i:i+1),fliplr(hostMedR.binEdges(i:i+1)));
    patch(xp,yp,cm(i,:),'facealpha',1,'edgecolor','k','edgealpha',1)
    hold on
    text(xt,yt,str,'fontsize',14,'interpreter','latex')
    %area(10.^(hostMed.binEdges(i:i+1)),10.^[yl(2) yl(2)],'facealpha',1);
end
 
h(1)=plot(xl(1)-1,yl(1)-1,'s','color',cmap(10,:),'markerfacecolor',cmap(10,:),...
'markersize',15,'DisplayName','$r_\mathrm{1/2,gas}$');
h(2)=plot(xl(1)-1,yl(1)-1,'s','color',cmap(3,:),'markerfacecolor',cmap(3,:),...
'markersize',15,'DisplayName','$r_\mathrm{star}$');

 hl=legend(h);
 set(hl,'Interpreter','latex','fontsize',16,'box','on','location','southEast');




% h(1)=plot(xdata,y1data,'o','markerfacecolor',cmdot(1,:),'color',cmdot(1,:),'markersize',4,'DisplayName','$r_\mathrm{star}$');
% h(2)=plot(xdata,y2data,'o','markerfacecolor',cmdot(2,:),'color',cmdot(2,:),'markersize',4,'DisplayName','$r_\mathrm{1/2,gas}$');
% 


box on

xlabelmine(mstarLab,18);
ylabelmine(ylab,18);
set(gca,'fontsize',14)

%% points and contours

ax1=axes;


data=dG+dS;
data(data==0)=-1;
data(dG==0 & dS==0)=-10;
imagesc(xl,yl,data)
caxis([-1.8 1.8])
colormap(cmap)
hold on


contour(popCont1.xx,popCont1.yy,popCont1.popContour,'ShowText','off','LineColor',concol,'linewidth',1,...
    'LevelList',lvl,'Fill','off','linestyle','-');

contour(popCont2.xx,popCont2.yy,popCont2.popContour,'ShowText','off','LineColor',concol,'linewidth',1,...
    'LevelList',lvl,'Fill','off','linestyle','-');

grid
set(gca,'Ydir','normal','Fontsize',14)


xlim(xl);ylim(yl);



set(ax1,'position',get(ax0,'position'));
linkaxes([ax0,ax1])
ax1.Visible = 'off';ax1.XTick = [];ax1.YTick = [];

%% rgas 

fname=['radiiPlot_' simDisplayName];
printout_fig(gcf,fname,'v','fig');

