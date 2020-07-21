%% plot 3-d mask 

siz=size(mask);


for i=1:siz(3)
    [r,c]=find(mask(:,:,i));
   
    z=i*ones(size(r));

    if(i==1)
        r1=r;
        c1=c;
        z1=z;
    else
        r1=cat(1,r1,r);
        c1=cat(1,c1,c);
        z1=cat(1,z1,z);
    end
end
    
    
    

