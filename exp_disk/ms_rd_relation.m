%% plot stripping catalog 

%load('mockCat_10000_xi.mat')

 hh=fspecial('gaussian',10,3);
%hh=fspecial('average',10);



[bir, binsize, xxlim,yylim]= histogram2d(log10(cata.Ms(maskCat)),log10(cata.rd(maskCat)),ones(size(cata.Ms(maskCat))),'xxlim',[8.8 12],'yylim',[-0.5 2.5],'len',[200 200]);
bird=bir(:,:,1)./length(cata.Ms(maskCat)).*100;
bs=imfilter(bird,hh);



figure
bb=zeros(size(bird));

for i=1:size(bb,1)
    for j=1:size(bb,2)
        mask=bs>=bs(i,j);
        bb(i,j)=sum(sum(bs(mask)));
    end
end

 xx=xxlim(1)+0.5*binsize(1):binsize(1):xxlim(2)-0.5*binsize(1);
 yy=yylim(1)+0.5*binsize(2):binsize(2):yylim(2)-0.5*binsize(2);
 lc=0:10:100;%lc(end+1)=99.9;
 bb(bb==max(max(bb)))=100;

 md=linspace(xxlim(1),xxlim(2),201);
rd=linspace(yylim(1),yylim(2),201);
 [mdMat,rdMat]=meshgrid(10.^md,10.^rd);
sigMat=mdMat./(2*pi.*rdMat.^2);
%load('MyColormaps','pop_hsv');
%colormap(pop_hsv)
 [C,h]= contourf(xx,yy,bb,lc);
 bar=colorbar;
 %clabel(C,h)
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
 hold on
[C1,hc1]=contour(md,rd,log10(sigMat),3:1:12,'-k','linewidth',2);
set(hc1,'ShowText','on','TextStep',get(hc1,'LevelStep')*2)

map=brewermap(256,'*YlOrRd');
map(257,:)=[1 1 1];
colormap(map);
%clabel(C1,hc1,'fontsize',18,'rotation',0)
%clabel(C1,'fontsize',14,'rotation',0)
xlim([8.8 12])
ylim([-0.5 2])
set(bar,'Fontsize',14,'Ytick',0:10:100)
set(get(bar,'Title'),'String','Population $[\%]$','Fontsize',14,'Interpreter','latex');
set(gca,'Ydir','normal','Fontsize',14)
%ylabelmine('$ \log(\Sigma_s) \, [\mathrm{M_\odot/kpc^2}]$',14)
ylabelmine('$ \log(R_{d}) \, [\mathrm{kpc}]$',14)
xlabelmine('$\log(M_s)\,[\mathrm{M_\odot}]$',14)

%% set up fitting stuff
 zz=log10(cata.rd);
 mm=log10(cata.Ms);
 indx=ceil((mm-xxlim(1))./(diff(xxlim)).*size(bs,2));
 indy=ceil((zz-yylim(1))./(diff(yylim)).*size(bs,1));
 
 wt=ones(size(mm));
 for i=1:length(cata.Ms)
     wt(i)=bs(indy(i),indx(i)).^2;
 end
 wt=wt./sum(wt);

 mas=(mm<10.5);
 mm1=mm(mas);
 zz1=zz(mas);
 wt1=wt(mas);
 wt1=wt1./sum(wt1);
 
 mm2=mm(~mas);
 zz2=zz(~mas);
 wt2=wt(~mas);
 wt2=wt2./sum(wt2);