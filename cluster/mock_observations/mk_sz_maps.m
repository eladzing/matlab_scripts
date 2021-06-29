%% generate SZ maps

%% perliminaries
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
units;
for lm=1:length(list)
    cl=list(lm);
    new_env(cl,'a1');
    
    %% generate projections for the boxes
    global NCELL;
    global hub;
    global CLUSTER;
    global zred;
    
    mask=true(NCELL,NCELL,NCELL);
    i1=floor(NCELL/4)+1;i2=i1+floor(NCELL/2)-1;
    
    unitFac=kb.*sigmaThomson/(me*cspeed^2);
    
    for i=1:4
        boxx=2^(i-1);
        cl=boxx./hub/NCELL.*Mpc;
        
        cub=RHOG(boxx).*(Ms/Mpc^3/mm).*T(boxx).*mask;
        
        szProj(i).xy=squeeze(cl.*trapz(cub,3)).*unitFac;
        szProj(i).yz=squeeze(cl.*trapz(cub,1)).*unitFac;
        szProj(i).xz=squeeze(cl.*trapz(cub,2)).*unitFac;
        szProj(i).cl=cl;
        szProj(i).boxx=boxx;
        if i==1
            mask(i1:i2,i1:i2,i1:i2)=false;
        end
    end
    
    clear cub mask
    %% generate integrated map
    len=szProj(end).boxx./szProj(1).boxx.*NCELL;
    
    projXY=zeros(len,len);
    projYZ=zeros(len,len);
    projXZ=zeros(len,len);
    
    indShift=(len-[1 2 4 8].*NCELL)/2;
    
    
    for k=1:length(szProj)
        
        boxx=szProj(k).boxx;
        for i=1:NCELL
            
            indx=indShift(k)+(((i-1)*boxx+1):boxx*i);
            for j=1:NCELL
                
                indy=indShift(k)+(((j-1)*boxx+1):boxx*j);
                
                projXY(indx,indy)=projXY(indx,indy) +  szProj(k).xy(i,j);
                projYZ(indx,indy)=projYZ(indx,indy) +  szProj(k).yz(i,j);
                projXZ(indx,indy)=projXZ(indx,indy) +  szProj(k).xz(i,j);
                
            end
        end
    end
    
    SZmap(lm).cluster=CLUSTER;
    SZmap(lm).zred=zred;
    SZmap(lm).projXY=projXY;
    SZmap(lm).projYZ=projYZ;
    SZmap(lm).projXZ=projXZ;
    
    
end
