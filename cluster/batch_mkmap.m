list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
mask=[ 1  1   1   1   1   1   1   1  1  1  1  1  1  1  1  1];


for i=1:length(list);
    if~mask(i) 
        continue;
    end
    
   new_env(sprintf('CL%d',list(i)),'csf','a1');
   
   for boxx=[8];  
%          circs=edge(i,:);
%          mkmap('Entropy',boxx,[-3 3],[1 1 1],1/256*6,'novel',-7,circs,[],'Edge by $\nabla S$','print','ent_a1_testedge1',printoutdir)  
%          mkmap('Temperature',boxx,[4 8],[1 1 1],1/256*6,'novel',-7,circs,[],'Edge by $\nabla S$','nprint','tmp_a1_testedge1',printoutdir)  
%          circs=edge2(i,:);
%          mkmap('Entropy',boxx,[-3 3],[1 1 1],1/256*6,'novel',-7,circs,[],'Edge by $max(S)$','print','ent_a1_testedge2',printoutdir)  
%          mkmap('Temperature',boxx,[4 8],[1 1 1],1/256*6,'novel',-7,circs,[],'Edge by $max(S)$','print','tmp_a1_testedge2',printoutdir)  
           
       mkmap('rps',8,[15 20],[1 1 1],boxx/256*6,'no',-4,[],1,'$z=0$','print','rps_a1')
   end
   
   close all
end