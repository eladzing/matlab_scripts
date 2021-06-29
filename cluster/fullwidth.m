tm=T(8);
ro=RHOG(8);
ltm=log10(tm);
lro=log10(ro);
logi=lro>11.5; % & ltm<6;
dm=2;
ttmp=rolg.*(tm.*ro);
rotmp=rolg.*ro;
imagesc(squeeze( mean(rotmp,dm)));
%imagesc(squeeze( sum(ttmp,dm)./sum(rotmp,dm)))
