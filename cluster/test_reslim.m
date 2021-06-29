list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];

res1=zeros(size(list));
res6=zeros(size(list));

for i=1:length(list)
    
    %boxx=1;
    new_env(sprintf('CL%d',list(i)),'csf','a1');
    rv=get_rvir;
    res1(i)=find_inner_reslim(1.0,5)./rv;
    
    new_env(sprintf('CL%d',list(i)),'csf','a06');
    rv=get_rvir;
    res6(i)=find_inner_reslim(1.0,5)./rv;
          
    
    
end

    
    
