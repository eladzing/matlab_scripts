
list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
mask=[ 1   0   1   0   1   1   0   1  1  1  0  1  1  1  1  1];

pos=zeros(length(list),4,6);
neg=zeros(length(list),4,6);
nc=zeros(length(list),4,6);

for i=1:16
   if~mask(i) 
        continue;
   end
    
   new_env(sprintf('CL%d',list(i)),'csf','a1');
   
   global CLUSTER;
    
   load(sprintf('mat_files/qdiss_np_%s_%s.mat','a1',CLUSTER));
   
   neg(i,:,:)=qdiss_array_np{2};
   nc(i,:,:)=qdiss_array_np{3};
end
   
 
   