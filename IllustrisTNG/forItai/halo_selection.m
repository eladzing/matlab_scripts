%% find the nearest neighnbor in units of ones own R200 
global illUnits

m200c=fofs.Group_M_Crit200.*illUnits.massUnit;

for i=1:fofs.count
    
    %tic
    newPos=double(illustris.utils.centerObject(fofs.GroupPos,fofs.GroupPos(:,i)));
    dist=sqrt(sum(newPos.^2,1));
    [mn,ix]=mink(dist,2);
    nearNeib.distance(i)=mn(2).*illUnits.lengthUnit;
    nearNeib.distanceNorm(i)=mn(2)./double(fofs.Group_R_Crit200(i));
    nearNeib.ID(i)=ix(2)-1;
    nearNeib.m200c(i)=m200c(ix(2));
    %toc
end

    
    
    