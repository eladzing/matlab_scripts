eta=logspace(-1,1,200);
zt=logspace(-1,0.7,200);
[rr,zz]=meshgrid(eta,zt);

gz=exp(-rr).*tanh(zz./2);
gr=0.5.*rr.*bfunc(rr);

gRatio=1-gr./(sqrt(gz.^2+gr.^2));


figure 

imagesc(log10(eta),log10(zt),gRatio)
set(gca,'Ydir','Normal','Fontsize',14,'Ytick',-1:0.2:0.7);

xlabelmine('$\log\left(R/R_s\right)$');
ylabelmine('$\log\left(z/z_s\right)$');
bar=colorbar;
barTitle(bar,'$\delta g/g$')
set(bar,'Ytick',0:0.1:1)
caxis([0 1]);
map=brewermap(256,'*RdYlBu');
colormap(map)