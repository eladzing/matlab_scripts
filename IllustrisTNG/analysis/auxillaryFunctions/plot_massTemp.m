function plot_massTemp(option,xx,yy,cc,treeSt,barLab,cmap,xl,yl,gasField,cax)

global tvirMean
global popCont

if xl(1)>=9
    mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
else
    mstarLab='$\log M_\mathrm{BH}\,[\mathrm{M_\odot}]$';
end


tempLab=sprintf('$\\log T_{\\mathrm{%s}} \\,[\\mathrm{K}]$',gasField);

%yl=[4 8];
%xl=[9 12.5];

if ~exist('cax','var')
    cax=[min(cc) max(cc)];
end


figure
set(gcf,'position',[1432 421 1000 750],'Color','w')

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
box on
hold on
errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',2)
% errorbar(tvirMean.xMedian,tvirMean.yMedian,...
%     tvirMean.yMedian-tvirMean.yQuarts(1,:),...
%     tvirMean.yQuarts(1,:)-tvirMean.yMedian,...
%     '--r','linewidth',1)
hb=colorbar;barTitle(hb,barLab,'fontsize',18)
set(hb,'fontsize',18);
colormap(cmap);caxis(cax);
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab,20);
ylabelmine(tempLab,20);
set(gca,'fontsize',18)

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',[10:10:100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];



end