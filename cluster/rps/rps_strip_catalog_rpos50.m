%% plot stripping catalog 

load('mockCat_1000.mat')

%etap=0.2;
mc=1e15;
fc=0.15;

etap=zeros(size(cat.sigma));


pval=rps_factor_expdisk('sigma_s',cat.sigma,'fd',cat.fg,...
    'fc',fc,'mc',mc,'etap',1);

eta=0.1:0.01:10;
%mstrip=zeros(size(cat.sigma));
pad=7;
tic
for i=1:length(etap)
    
f1=disk_force_reduced(eta,'beta',cat.beta(i),'fg',cat.fg(i),...
    'BT',cat.BT(i));
[f1Max, id]=max(f1);

f1(1:id)=f1Max;

mf1=exp_disk_mass(eta,cat.beta(i));

fi=interp1(mf1,f1,0.5);

etap(i)=(fi./pval(i)).^(-0.5);

        
    
end
toc

%mstrip=(1-mstrip).*100;

% figure
% subplot(3,1,1:2)
% plot(log10(cat.sigma),mstrip,'.')
% xlim([6 11])
% grid
% ylabelmine('$M_{strip}\,[\%]$',14)
% set(gca,'Fontsize',14,'box','on');
% 
% subplot(3,1,3)
% [nh xo]=hist(log10(cat.sigma),20);
% bar(xo,nh./length(mstrip).*100);
% set(gca,'Fontsize',14,'box','on');
% xlim([6 11])
% xlabelmine('$\Sigma_s\,[\mathrm{M_\odot/kpc^{2}}]$',14)
% ylabelmine('$\%$',14)





 hh=fspecial('gaussian',6,3);



[bir, binsize, xxlim,yylim]= histogram2d(log10(cat.sigma),log(etap),ones(size(etap)),'xxlim',[6 11],'yylim',[-3 2]);
bird=bir(:,:,1)./length(etap).*100;
bss=imfilter(bird,hh);
%bs=imfilter(bird,hh,'full');
 bs=bss./sum(sum(bss)).*100;

%figure
bb=zeros(size(bird));

for i=1:size(bb,1)
    for j=1:size(bb,2)
        mask=bs>=bs(i,j);
        bb(i,j)=sum(sum(bs(mask)));
    end
end
% contourf(bb)
% contourf(bb,10:10:100)
% contourf(bb,1:10:100)
 xx=xxlim(1)+0.5*binsize(1):binsize(1):xxlim(2)-0.5*binsize(1);
 yy=yylim(1)+0.5*binsize(2):binsize(2):yylim(2)-0.5*binsize(2);
 lc=[0 10 20 30 40 50 60 70 80 90 100] ;%    0:10:100;%lc(end+1)=99.9;
 bb(bb==max(max(bb)))=100;

 createfigure_stripped_catalog_with_cumsum(xx, yy, bb, cumsum(sum(bird,2)))
 
% [C,h]= contourf(xx,yy,bb,lc);
% clabel(C,h)
% bar=colorbar;
% set(bar','Fontsize',14)
% set(get(bar,'Title'),'String','Population $[\%]$','Fontsize',14,'Interpreter','latex');
% set(gca,'Ydir','normal','Fontsize',14)
% xlabelmine('$ \log(\Sigma_s) \, [\mathrm{M_\odot/kpc^2}]$',14)
% ylabelmine('Stripped Mass $[\%]$',14)
