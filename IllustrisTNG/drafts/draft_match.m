matched=zeros(size(qGroupInd));

massFac=0.20;
for j=1:100
    tempList=oGroupInd;
    
    for i=1:length(qGroupInd)
        range=galMass(qGroupInd(i)).*(1+[-1 1].*massFac);
        
        mm=galMass(tempList)>=range(1) & galMass(tempList)<=range(2);
        
        if sum(mm)>1
            
            cand=tempList(mm);
            
            
            ran=rand(size(cand));
            
            [~,choice]=max(ran);
            
            matched(i)=cand(choice);
            
            tempList=tempList(tempList~=cand(choice));
            
               
            
            
        end
    end
    
    maxfac2(j)=max(abs(galMass(matched)./galMass(qGroupInd)-1));
    meanfac2(j)=mean(abs(galMass(matched)./galMass(qGroupInd)-1));
    medfac2(j)=median(abs(galMass(matched)./galMass(qGroupInd)-1));
end
%%
figure;loglog(galMass(qGroupInd),galMass(matched),'.')

figure;semilogx(galMass(qGroupInd),galMass(matched)./galMass(qGroupInd),'.')