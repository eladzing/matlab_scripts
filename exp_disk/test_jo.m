


qMask=cata.lambdaMask & cata.qmin>=1;

mt=cata.Ms.*(1+cata.fb);

massMask=log10(mt)>=9 & log10(mt)<=9.5;

msk=qMask & massMask;

ssfr=ssfr_proxy_int(cata.Ms(msk),cata.rd(msk),cata.fg(msk),cata.beta(msk),cata.fb(msk));

re=[];
for i=1:length(cata.Ms)
    if ~msk(i)
        continue
    end
    re(end+1)=find_re(cata.fb(i),cata.xi(i)).*cata.rd(i);
end
re=re';
sig=0.5.*mt(msk)./(pi.*re.^2);

%sig=cata.sigmaeff(msk);


[bir, binsize, xxlim,yylim]= histogram2d(log10(sig),log10(ssfr),ones(size(sig)),'len',[50 50]);
bird=bir(:,:,1)./length(sig).*100;
hh=fspecial('gaussian',3,3);
bs=imfilter(bird,hh);
figure
bb=zeros(size(bird));

for i=1:size(bb,1)
    for j=1:size(bb,2)
        mask=bs>=bs(i,j);
        bb(i,j)=sum(sum(bs(mask)));
    end
end

map=brewermap(256,'*YlOrRd');
 xx=xxlim(1)+0.5*binsize(1):binsize(1):xxlim(2)-0.5*binsize(1);
 yy=yylim(1)+0.5*binsize(2):binsize(2):yylim(2)-0.5*binsize(2);
 lc=0:10:100;%lc(end+1)=99.9;
 bb(bb==max(max(bb)))=100;
 [C,h]= contourf(xx,yy,bb,lc);
 colormap(map)
 
 bar=colorbar;
 barTitle(bar,'$\%$');
 set(gca,'fontsize',14,'box','on')
 
 xlabelmine('$\log(\Sigma_e)\,[\mathrm{M_\odot/kpc^2}$]')
 ylabelmine('SSFR $[\mathrm{yr^{-1}}]$')
 
 

