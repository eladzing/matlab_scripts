function plot_bhratEnt2(option,xx,yy,cc,treeSt,barLab,cmap,xl,yl,popCont,qCont,sfCont,gasField,cax,hostm)


%   Detailed explanation goes here        

%global tvirMean
%global popCont

global illUnits
global DRACOFLAG

if ~DRACOFLAG
    tagFont=16;
    axFont=18;
    bigTagFont=18;
    lw=3;
else
    tagFont=30;
    axFont=30;
    bigTagFont=34;
    lw=5;
end

switch gasField
    case 'Gal'
        gTag={'Galactic Gas' '$r<r_\mathrm{gal}$' '(w/o SF gas)'};
    case 'CGM'
        gTag={'Inner CGM' '$r_\mathrm{gal}<r<r_\mathrm{1/2,gas}$'};
    case 'Out'
        gTag={'Outer CGM' '$r>r_\mathrm{1/2,gas}$'};
    case 'Sub'
        gTag={'CGM' '$r>r_\mathrm{gal}$'};
end
 


bhLab='log Cum. kinetic/thermal energy injection ratio';
entLab='$\log K\,[\mathrm{keV\,cm^2}]$';




%yl=[-3 2.5];
%xl=[9 12.5];

if ~exist('cax','var')
    cax=[min(cc) max(cc)];
end


figure
set(gcf,'position',[1432 421 1000 750],'Color','w')

%% background with halo masses

ax0=axes; 
xlim(xl);ylim(yl);

len=length(hostm.binEdges)-1;
cm=brewermap(len+1,'RdPu');

yp=[yl(1) yl(1) yl(2) yl(2)];
yt=yl(2)-0.05*diff(yl);
for i=1:len
%     if i==1 
%     str=sprintf('$M_\\mathrm{H}=10^{%3.1f}$',hostm.yMedian(i));
%       %xt=xl(1)+0.1.*diff(hostm.binEdges(i:i+1));
%     else
%         
%         %xt=hostm.binEdges(i)+0.3.*diff(hostm.binEdges(i:i+1));
%        
%     end
    str=sprintf('$10^{%3.1f}$',hostm.yMedian(i));
    xt=hostm.binEdges(i)+0.3.*diff(hostm.binEdges(i:i+1));
    xp=cat(2,hostm.binEdges(i:i+1),fliplr(hostm.binEdges(i:i+1)));
    patch(xp,yp,cm(i,:),'facealpha',0.2,'edgecolor','k','edgealpha',0.4)
    hold on
    text(xt,yt,str,'fontsize',tagFont,'interpreter','latex')
    %area(10.^(hostm.binEdges(i:i+1)),10.^[yl(2) yl(2)],'facealpha',1);
end
grid
box on

% gas tag


% gas tag
xfac=0.08; yfac=0.95;
text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),gTag{1},...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
    'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold')

xfac=0.08; yfac=0.89;
text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),gTag{2},...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
    'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold')

if strcmp(gasField,'Gal')
    xfac=0.24; yfac=0.88;
    text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),gTag{3},...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
        'Interpreter','latex','fontsize',tagFont-4)
end


% redshift tag
if illUnits.zred ==0
    zTag=sprintf('z=%i',illUnits.zred);
else
    zTag=sprintf('z=%3.2f',illUnits.zred);
end

xfac=0.9; yfac=0.05;
text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),zTag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
    'Interpreter','latex','fontsize',tagFont,'fontweight','bold')


xlabelmine(bhLab,tagFont);
ylabelmine(entLab,tagFont);
set(gca,'fontsize',axFont)
%ax0.Visible = 'off';ax0.XTick = [];ax0.YTick = [];

%% plot color histogram

ax1=axes;
switch(lower(option))
    case {'scatter','scat','points'}
        
        scatter(xx,yy,18,cc)
    case {'tree','cells'}
        celVal=points2tree(cc,treeSt,'median');
        plot2dTreeMap(celVal,treeSt,'cmap',cmap,'fig',gcf,'minmax',cax);
    otherwise
        error('plot_massTemp - illegal option: %s',option)
end

hold on

contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',5:10:95,'Fill','off','linestyle','-');

qCol=[0.89,0.1,0.11];
[~,h(1)]=contour(qCont.xx,qCont.yy,qCont.popContour,'ShowText','off',...
    'LineColor',qCol,'linewidth',8,...
    'LevelList',90,'Fill','off','linestyle','--',...
    'DisplayName','Quenched');

sCol=[0.1,0.1,0.89];
[~,h(2)]=contour(sfCont.xx,sfCont.yy,sfCont.popContour,'ShowText','off',...
    'LineColor',sCol,'linewidth',8,...
    'LevelList',90,'Fill','off','linestyle','--',...
    'DisplayName','Star-forming');


xfac=0.8; yfac=0.95;
text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),'Quenched',...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
    'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold','color',qCol)

xfac=0.33; yfac=0.045;
text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),'Star-forming',...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
    'Interpreter','latex','fontsize',bigTagFont,'fontweight','bold','color',sCol)


if iscell(barLab)
    barCell=barLab;
    if length(barCell)==1
        barCell{2}='out';
    end
else
    barCell{1}=barLab;
    barCell{2}='out';
end

hb=colorbar;
switch(lower(barCell{2}))
    case 'in'
        
        %hb=colorbar;%barTitle(hb,barLab,'fontsize',tagFont)
        set(hb,'fontsize',axFont,'axisLocation','in',...
            'position',[0.88 0.13 0.02 0.58]);
        
        xfac=1.025; yfac=0.5;
        text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),barCell{1},...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
            'Interpreter','latex','fontsize',tagFont,'color','k','rotation',270)
        
    case 'out'
        
        %barTitle(hb,barLab,'fontsize',tagFont)
        set(hb,'fontsize',axFont);
        xfac=0.82; yfac=1.05;
        text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),barCell{1},...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
            'Interpreter','latex','fontsize',tagFont,'color','k','rotation',0)
        
    case 'innone'
        set(hb,'fontsize',axFont,'axisLocation','in',...
            'position',[0.88 0.13 0.02 0.58]);
    case 'none'
        
    otherwise
        error('PLOT_MASSENT3 - illegal bar option: %s',barCell{2})
end

colormap(cmap);caxis(cax);
xlim(xl);ylim(yl);
set(gca,'fontsize',axFont)
set(ax1,'position',get(ax0,'position'));

linkaxes([ax0,ax1]) 
ax1.Visible = 'off';ax1.XTick = [];ax1.YTick = [];

% 
% ax2=axes;
% contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','off','LineColor',[0 0 0],...
%     'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
% xlim(xl);ylim(yl);
% set(ax2,'position',get(ax1,'position'));
% linkaxes([ax0,ax1,ax2])
% ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
% ax1.Visible = 'off';ax1.XTick = [];ax1.YTick = [];

end

