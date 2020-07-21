
splashBack1=false(size(centralHist));
splashBack2=splashBack1;

dThresh=0.05;

for i=1:length(done)
    
    if done(i)

        %% not central sub condition 
        satInd=find(centralHist(i).isCentral==0);
        
        cond1=length(satInd)>1;
        cond2=any(diff(satInd)==1);
        
        splashBack1(i)=cond1 & cond2;
        
        
        %% distance from center cond
        
        rpos=centralHist(i).radiusToHost./centralHist(i).hostR200;
        satInd2=find(rpos>dThresh);
        
        cond1=length(satInd2)>1;
        cond2=any(diff(satInd2)==1);
        
        splashBack2(i)=cond1 & cond2;
    end
    
end

        
        