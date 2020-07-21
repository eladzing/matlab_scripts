lm=8:0.01:9;
ms=10.^lm;
mdList=[0.05 0.1 0.15];
mdList=mdList(2);
sampleSize=1;
rd=zeros(size(ms));
rdOld=rd;
%wtbarLen=sampleSize;



lamdaArr=lambda_prime(ms);

tic
res=rscale_mmw_array(ms,'md',0.1,'lambda',lamdaArr,'noshow');
rd=res.rd;
toc

tic 
for i=1:length(ms)
    res=rscale_mmw(ms(i),'md',0.1,'lambda',lamdaArr(i),'noshow');
    rdOld(i)=res.rd;
end
toc 
