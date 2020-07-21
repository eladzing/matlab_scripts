%% generate multi scale pot maps

%% perliminaries
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
list=list(10);
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
    
    unitFac=1; %kb.*sigmaThomson/(me*cspeed^2);
    
    for i=1:4
        boxx=2^(i-1);
        cel=boxx./hub/NCELL.*1000; % inkpc 
        [pdm,pgs,~,~]=poteng(boxx);
        cub=((pdm+pgs)./(1000.*cel).^3).*mask; % in units of (km/sec)^2/kpc^3
        
        potProj(i).xy=squeeze(cel.*trapz(cub,3)).*unitFac;% in units of (km/sec)^2/kpc^2
        potProj(i).yz=squeeze(cel.*trapz(cub,1)).*unitFac;
        potProj(i).xz=squeeze(cel.*trapz(cub,2)).*unitFac;
        potProj(i).cel=cel;
        potProj(i).boxx=boxx;
         if i==1
             mask(i1:i2,i1:i2,i1:i2)=false;
         end
    end
    
    clear cub mask
    %% generate integrated map
    len=potProj(end).boxx./potProj(1).boxx.*NCELL;
    
    projXY=zeros(len,len);
    projYZ=zeros(len,len);
    projXZ=zeros(len,len);
    
    indShift=(len-[1 2 4 8].*NCELL)/2;
    
    
    for k=1:length(potProj)
        
        boxx=potProj(k).boxx;
        for i=1:NCELL
            
            indx=indShift(k)+(((i-1)*boxx+1):boxx*i);
            for j=1:NCELL
                
                indy=indShift(k)+(((j-1)*boxx+1):boxx*j);
                
                projXY(indx,indy)=projXY(indx,indy) +  potProj(k).xy(i,j);
                projYZ(indx,indy)=projYZ(indx,indy) +  potProj(k).yz(i,j);
                projXZ(indx,indy)=projXZ(indx,indy) +  potProj(k).xz(i,j);
                
            end
        end
    end
    
    potMap(lm).cluster=CLUSTER;
    potMap(lm).zred=zred;
    potMap(lm).projXY=projXY;
    potMap(lm).projYZ=projYZ;
    potMap(lm).projXZ=projXZ;
    
    
end
