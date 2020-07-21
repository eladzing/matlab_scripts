bp=illustris.set_env(100);

global massUnit
global LBox

snaps=[99 84 72 63];

for i=1:length(snaps)
    
    snp=snaps(i);
    fprintf(' *** beginning snap %i ***\n',snp);
    
    if readFlag
        subs=illustris.groupcat.loadSubhalos(bp,snp);
        fofs=illustris.groupcat.loadHalos(bp,snp);
    end
    
    fprintf(' *** finished reading ***\n');
    
    
    %% 1. Host’s Group_M_Crit200 >= 1E13 Msun
    maskHostMass=(fofs.Group_M_Crit200(subs.SubhaloGrNr+1).*massUnit)>1e13;
    
    %% SubhaloMassinRadStar >= 1E9.5 Msun
    maskGalMass=(subs.SubhaloMassInRadType(illustris.partTypeNum('star')+1,:).*massUnit)>10^9.5;
    
    %% SubhaloMassInRadDM > 0.1 * SubhaloMassInRad
    maskDmFrac=subs.SubhaloMassInRadType(illustris.partTypeNum('dm')+1,:)./subs.SubhaloMassInRad>0.1;
    %maskDmFrac=subs.SubhaloMassType(illustris.partTypeNum('dm')+1,:)./subs.SubhaloMass  >0.1;
    
    %% SubhaloMassInRadGas > 0
    maskGasMass=subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,:)>0;
    
    %% central condition 
    subIndx=0:subs.count-1;
    isCentral=subIndx==fofs.GroupFirstSub(subs.SubhaloGrNr+1);
    
    
    %% dist(SubhaloPos, GroupPos) <= Group_M_Crit200 or dist(Subhalo, Host) <= 2xGroup_M_Crit200 at z = z_select
    
    dist=zeros(1,subs.count);
    for j=1:3
        dx=abs(subs.SubhaloPos(j,:)-fofs.GroupPos(j,subs.SubhaloGrNr+1));
        m1=dx>0.5.*LBox;
        dx(m1)=dx(m1)-LBox;
        
        dist=dist+dx.^2;
    end
    dist=sqrt(dist)./fofs.Group_R_Crit200(subs.SubhaloGrNr+1);
    
    maskDist1=dist<1;
    
    maskDist2=dist<2;
    
    %% calculate numbers
    
    basePopulation1(i)=sum(maskHostMass & maskGalMass & maskDmFrac & maskGasMass & maskDist1 & ~isCentral);
    
    basePopulation2(i)=sum(maskHostMass & maskGalMass & maskDmFrac & maskGasMass & maskDist2 & ~isCentral);
    
    
    
    
end
