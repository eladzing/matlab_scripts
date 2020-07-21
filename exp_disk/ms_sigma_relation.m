%% plot stripping catalog 

load('mockCat_1000.mat')

 hh=fspecial('gaussian',20,5);
%hh=fspecial('average',10);

[bir, binsize, xxlim,yylim]= histogram2d(log10(cat.Ms),log10(cat.sigma),ones(size(cat.Ms)),'xxlim',[8 12],'yylim',[6 11]);
bird=bir(:,:,1)./length(cat.Ms).*100;
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
[C,h]= contourf(xx,yy,bb,lc);
%clabel(C,h)
bar=colorbar;
set(bar','Fontsize',14)
set(get(bar,'Title'),'String','Population $[\%]$','Fontsize',14,'Interpreter','latex');
set(gca,'Ydir','normal','Fontsize',14)
ylabelmine('$ \log(\Sigma_s) \, [\mathrm{M_\odot/kpc^2}]$',14)
xlabelmine('$\log(M_s) $[\mathrm{M_\odot}]$',14)
