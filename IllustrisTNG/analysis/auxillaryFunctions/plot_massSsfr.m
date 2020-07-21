function hf=plot_massSsfr(option,xx,yy,cc,treeSt,barLab,cmap,cax)

%global tvirMean
%global popCont
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
ssfrLab='$\log \mathrm{sSFR}\,[\mathrm{yr^{-1}}] $';

xl=[9 12.5];
yl=[-16.5 -8];

if ~exist('cax','var')
    cax=[min(cc) max(cc)];
end

hf=figure;
set(gcf,'position',[1432 421 1000 750])

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

hb=colorbar;barTitle(hb,barLab,'fontsize',16)
set(hb,'fontsize',14);
colormap(cmap);caxis(cax);
xlim(xl);ylim(yl);
grid
xlabelmine(mstarLab,18);
ylabelmine(ssfrLab,18);
set(gca,'fontsize',14)

ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','of','LineColor',[0 0 0],...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];


end