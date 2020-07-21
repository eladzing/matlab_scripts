list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];

for i=1:length(list)
    
    %boxx=1;
    cl=sprintf('CL%d',list(i));
   
    [qb1,~,~,qb8] = cooling_profs(cl,'csf','a1',[0 0 0],0.02);
    
       
    rv=get_rvir;
    figure    
    r1=(1:256).*1./0.7./256./2;
    r8=(1:256).*8./0.7./256./2;
    loglog(r1./rv,qb1,r8./rv,qb8)
    xlim([0.01 0.1])
    
    
end

    
