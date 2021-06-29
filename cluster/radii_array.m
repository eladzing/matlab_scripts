list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

rr=zeros(length(list),2);

for i=1:length(list)
    
    new_env(list(i));
    
    rr(i,1)=get_rvir(500);
    rr(i,2)=get_rvir(200);
end