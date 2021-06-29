%% make bird 

new_env(6)
g1=RHOG(1);g2=RHOG(2);g4=RHOG(4);g8=RHOG(8);
t1=T(1);t2=T(2);t4=T(4);t8=T(8);
mask=true(size(g1));
mask(65:192,65:192,65:192)=false;
mask0=true(size(g1));
vol=([1 2 4 8]/0.7/256).^3;

m1=g1(mask0).*vol(1);
m2=g2(mask).*vol(2);
m4=g4(mask).*vol(3);
m8=g8(mask).*vol(4);

gg2=g2(mask);
gg1=g1(mask0);
gg4=g4(mask);
gg8=g8(mask);
tt1=t1(mask0);
tt2=t2(mask);
tt4=t4(mask);
tt8=t8(mask);
gg=cat(1,gg1,gg2,gg4,gg8);
tt=cat(1,tt1,tt2,tt4,tt8);
mm=cat(1,m1,m2,m4,m8);
[bird, binsize, xxlim,yylim]= histogram2d(log10(tt),log10(gg),mm);

cmap=brewermap(256,'YlOrRd');
imagesc(xxlim,yylim,log10(squeeze(bird(:,:,1))))
set(gca,'Ydir','Normal');
colormap(cmap);
colorbar
xlabelmine('$\log(T)\,[\mathrm{K}]$')
ylabelmine('$\log(\rho)\,[\mathrm{M_\odot\, Mpc^{-3}}]$')

