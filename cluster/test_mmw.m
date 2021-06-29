lm=8:0.05:11;
ms=10.^lm;
mdList=[0.05 0.1 0.15];
mdList=mdList(2);
sampleSize=10;
rd=zeros(length(ms),length(mdList),sampleSize);
rdOld=rd;
wtbarLen=sampleSize;

tic
hwb=waitbar(0);
cnt=0;
for i=1:sampleSize  %loop over sample
    %lamdaArr=zeros(size(ms));
    lamdaArr=lambda_prime(ms);
    if mod(i,10)==0
        waitbar(i/wtbarLen,hwb)
    end
    for j=1:length(mdList)
        md=mdList(j);
        %tic
        %res=rscale_mmw_array(ms,'md',md,'lambda',lamdaArr);
        %rd(:,j,i)=res.rd;
        %toc
        
        
        for k=1:length(ms);
            res=rscale_mmw(ms(k),'md',md,'lambda',lamdaArr(k),'noshow');
        %    if ~res.converged
        %        error('not convereged')
        %    end
            rdOld(k,j,i)=res.rd;
            
            %pause
        end
       
    end
end
 toc
close(hwb)



