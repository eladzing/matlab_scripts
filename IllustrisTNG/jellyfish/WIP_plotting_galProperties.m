
propSt 




h=[];
if ~exist('cax','var')
    cax=[min(cc) max(cc)];
end



%% 

xx=log10(galProps.galStellarMass);
yy=log10(galProps.hostM200c);

xl0=[min(xx(mask100)) max(xx(mask100))].*[0.99 1.01];
yl0=[min(yy(mask100)) max(yy(mask100))].*[0.99 1.01];
%% prepare data for plotting  
[bird, binsize, xl,yl]= histogram2d(xx(mask100 & ~maskJF),yy(mask100 & ~maskJF),...
ones(size(xx(mask100))),'xlim',xl0,'ylim',yl0,'len',50);

filt=fspecial('disk',15);
popCont=plot_population_contour(xx(maskJF & mask100),yy(maskJF & mask100),'smooth',filt,'noplot');

jfscore=ones( size(xx(mask100)));
jfscore(~maskJF)=0;


%% identify points beyond contour
xxx=xx(mask100 & maskJF);
yyy=yy(mask100 & maskJF);

lx=max(diff(popCont.xx));ly=max(diff(popCont.yy));
x0=popCont.xx(1)-0.5.*lx; y0=popCont.yy(1)-0.5.*ly;

indx=ceil((xxx-x0)./lx); indy=ceil((yyy-y0)./ly);
for i=1:length(xxx)
    score(i)=popCont.popContour(indy(i),indx(i));
end


cmap=brewermap(256,'Reds');
cmap=inferno(256);
qCol=brewermap(8,'Set1');

jfCol=[0.1,0.1,0.89];
njfCol=[0.89,0.1,0.11];
cols=cat(1,njfCol,jfCol);
%% background with full sample grayscale
figure
figure('position',[1432 421 1000 750],'Color','w')

%% underlying hist 

hh=scatterhist(xx(mask100),yy(mask100),'Group',jfscore(mask100),...
    'location','northeast','direction','out','legend','off','color',cols,...
    'markersize',1);

%0.6950    0.6960
set(hh(1),'position',[0.1000    0.1000   0.7 0.7 ],'fontsize',14);
hh(1).YLabel.String='';
hh(1).XLabel.String='';
    %'XTick',[],'YTick',[],'YAxisLocation','left','XAxisLocation','bottom');
set(hh(2),'position',[0.80 0.1 0.119806451612903 0.75]);
set(hh(3),'position',[0.1 0.80 0.75 0.132374233128834]);
hh(2).Children(1).LineWidth=1.5;
hh(2).Children(2).LineWidth=1.5;
hh(3).Children(1).LineWidth=1.5;
hh(3).Children(2).LineWidth=1.5;
xlim(xl);ylim(yl);

%% map & contout 
ax1=axes;
pos=get(hh(1),'position');
set(ax1,'position',pos);

imagesc(xl,yl,squeeze(bird(:,:,1)))
set(gca,'ydir','normal','fontsize',14)

colormap(cmap)

%% JF


hold on
%qCol=[0.89,0.1,0.11];
sCol=[0.1,0.1,0.89];
[~,h(1)]=contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','off','LineColor',[0 0 1],...
 'LineColor',sCol,'linewidth',2,...
    'LevelList',[98 75:-25:5],'Fill','off','linestyle','-',...
    'DisplayName','Jellyfish');    


plot(xxx(score>98),yyy(score>98),'^','color',sCol,'markersize',6.5,'linewidth',1.5,'MarkerFaceColor',qCol(2,:),...
    'Displayname','JF beyond 98 percentile');

%legend(h)

xfac=0.83; yfac=0.08;
text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),'TNG100',...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
    'Interpreter','latex','fontsize',17,'fontweight','bold','color','k')



xlim(xl);ylim(yl);

xlabelmine('log Stellar Mass',16)
ylabelmine('log host $M_\mathrm{200}$',16)


linkaxes([hh(1),ax1])  %linkaxes([ax0,ax1,ax2])
%ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
%ax1.Visible = 'off';ax1.XTick = [];ax1.YTick = [];
%hh(1).Visible = 'off';hh(1).XTick = [];hh(1).YTick = [];
