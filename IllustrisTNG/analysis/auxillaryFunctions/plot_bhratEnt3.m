function plot_bhratEnt2(option,xx,yy,cc,treeSt,barLab,cmap,xl,yl,popCont,gasField,cax,hostm)
%   Detailed explanation goes here

%global tvirMean
%global popCont

switch gasField
    case 'Gal'
        tag='ISM Gas';
    case 'CGM'
        tag='Inner CGM';
    case 'Out'
        tag='Outer CGM';
    case 'Sub'
        tag='Sub-halo Gas';
end
        
 


bhLab='$\log E_\mathrm{LAM}/E_\mathrm{HAM}$';
entLab='$\log S\,[\mathrm{KeV\,cm^2}]$';




%yl=[-3 2.5];
%xl=[9 12.5];

if ~exist('cax','var')
    cax=[min(cc) max(cc)];
end


figure
set(gcf,'position',[1432 421 1000 750],'Color','w')

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
    text(xt,yt,str,'fontsize',16,'interpreter','latex')
    %area(10.^(hostm.binEdges(i:i+1)),10.^[yl(2) yl(2)],'facealpha',1);
end
grid
box on

xfac=0.2; yfac=0.9;
text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),tag,...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
'Interpreter','latex','fontsize',28,'fontweight','bold')


xlabelmine(bhLab,20);
ylabelmine(entLab,20);
set(gca,'fontsize',18)
%ax0.Visible = 'off';ax0.XTick = [];ax0.YTick = [];



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
% hold on
% errorbar(tvirMean.xMedian,tvirMean.yMean,tvirMean.yStanDev/2,'-k','linewidth',1)
hb=colorbar;barTitle(hb,barLab,'fontsize',18)
set(hb,'fontsize',20);
colormap(cmap);caxis(cax);
xlim(xl);ylim(yl);

set(ax1,'position',get(ax0,'position'));




ax2=axes;
contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','off','LineColor',[0 0 0],...
    'LevelList',[5:10:95 100],'Fill','off','linestyle','-');
xlim(xl);ylim(yl);
set(ax2,'position',get(ax1,'position'));
linkaxes([ax0,ax1,ax2])
ax2.Visible = 'off';ax2.XTick = [];ax2.YTick = [];
ax1.Visible = 'off';ax1.XTick = [];ax1.YTick = [];

end

