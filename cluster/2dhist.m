ltm=log10(t');
lro=log10(ro');
tro(:,2)=ltm(1:256^3);
tro(:,1)=lro(1:256^3);
n=hist3(tro,[100 100]);
imagesc(n)
