% plot all simulation cells on the P-ro plane and find polytropic index

new_env(7)

% box1

for k=[1 2 4 8]
    
    rog=RHOG(k);
    tmp=T(k);
    
    mask=true(size(rog));
    if k==1
        mask(124:133,124:133,124:133)=false;
    else
        mask(65:192,65:192,65:192)=false;
    end
    
    rog=rog(mask);
    tmp=tmp(mask);
    
    pre=rog.*tmp;
    mass=rog.*get_cellsize(k,'Mpc').^3;
    clear tmp mask
    
    if k==1
        roTot=rog;
        preTot=pre;
        massTot=mass;
    else
        roTot=cat(1,roTot,rog);
        preTot=cat(1,preTot,pre);
        massTot=cat(1,massTot,mass);
    end
    
end

[bird, binsize, xxlim,yylim]= histogram2d(log10(roTot),log10(preTot),massTot);
