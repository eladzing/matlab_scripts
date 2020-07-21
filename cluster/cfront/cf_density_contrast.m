% plot cold fron density contrast 

gm=5/3;


m1=logspace(0,1,500);   % 1:0.01:6;
m2=m1;

[m1Mat,m2Mat]=meshgrid(m1,m2);


gf=2/(gm-1);

qf=(1+gf./m2Mat.^2)./(1+gf./m1Mat.^2);

% contourf(m1,m2,log10(qf),-0.5:0.05:0.5)
% bar=colorbar;
% 
% xlabelmine('$\mathcal{M}_1$')
% ylabelmine('$\mathcal{M}_2$')

%load('MyColormaps','newJet');
map = brewermap(256,'*RdBu');

hf=figure;
imagesc(log10(m1),log10(m2),log10(qf))
%imagesc(m1,m2,log10(qf))
hold on
%cl=log10([ 0.33 0.5 1 1.5 2 3 ]);
cl=log10([0.33 0.5 0.9 1 1.1  2 3 ]);
[C,h]=contour(log10(m1),log10(m2),log10(qf),cl,'k');
%clabel(C,'manual','fontsize',12)
%clabel(C,'manual','fontsize',12)
lab={'\alpha=0.33','\alpha=0.5','\alpha=0.9','\alpha=1','\alpha=1.1','\alpha=2','\alpha=3'};
tc=clabel(C,'manual','fontsize',12);
for i=2:2:length(tc);
    set(findobj(tc,'String',get(tc(i),'String')),'String',lab(floor(i/2)))
end
%set(h,'ShowText','on','Text',{'1/3','1/2','1','1.1','1.5','2','3'})
caxis(log10([0.25 4]));
colormap(map);
set(gca,'Ydir','Normal','Fontsize',14)

bar=colorbar;
barTitle(bar,'$\log(\alpha)$')
set(bar,'fontsize',14);
xlabelmine('$\log(\mathcal{M}_1)$')
ylabelmine('$\log(\mathcal{M}_2)$')

%printout_fig(hf,'qfactor_cfront','subdir','cfront')
