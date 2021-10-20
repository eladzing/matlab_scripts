
propSt 




h=[];
if ~exist('cax','var')
    cax=[min(cc) max(cc)];
end

figure
set(gcf,'position',[1432 421 1000 750],'Color','w')


%% 

xx=log10(galProps.galStellarMass);
yy=log10(galProps.hostM200c);

[bird, binsize, xl,yl]= histogram2d(xx(mask100 & ~maskJF),yy(mask100 & ~maskJF),ones(size(xx(mask100))),'len',50);

filt=fspecial('disk',15);
popCont=plot_population_contour(xx(maskJF & mask100),yy(maskJF & mask100),'smooth',filt,'noplot');

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
qCol=brewermap(8,'Set1');
%% background with full sample grayscale
figure

%ax0=axes;
xlim(xl);ylim(yl);


imagesc(xl,yl,squeeze(bird(:,:,1)))
set(gca,'ydir','normal')

colormap(cmap)

%% JF

%ax1=axes;
hold on
%qCol=[0.89,0.1,0.11];
sCol=[0.1,0.1,0.89];
[~,h(1)]=contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','off','LineColor',[0 0 1],...
 'LineColor',sCol,'linewidth',1.8,...
    'LevelList',[98 75:-25:5],'Fill','off','linestyle','-',...
    'DisplayName','Jellyfish');    


plot(xxx(score>98),yyy(score>98),'^','color',sCol,'markersize',6.5,'linewidth',1.5,'MarkerFaceColor',qCol(2,:),...
    'Displayname','JF beyond 98 percentile');

legend(h)

xlim(xl);ylim(yl);
xlabelmine