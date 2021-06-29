load(sprintf('mat_files/cooling_profiles_%s.mat','a1'));

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

ll=length(list);
outpath='/home/eladzing/work/clusters/printout';
for i=1:ll 
    clname=sprintf('CL%d',list(i))
    
    qb=coolprofs{i,6};
    rq=qb(1,:);
    q=qb(2,:);
    clear qb;
    [inedg rinedg]=find_inner_cooling_edge(rq,q,clname);
    
    coolprofs{i,7}=[inedg rinedg];
    
    saveas(gcf,sprintf('%s/%s_coolprof_testderiv.png',outpath,clname));
    close gcf
end
%     
%     
%     
%     figure;
%     loglog(r1,q1,r2,q2./2,r4,q4./4,r8,q8./8)
%     xlabel('r [Mpc]')    
%     ylabel('q [M_\odot*(km/sec)^2]')
%     title(sprintf('%s Cooling Profiles normalized',clname)); 
%     saveas(gcf,sprintf('%s/%s_coolprof2.png',outpath,clname));
%     
% %     figure;
% %     loglog(rq,cumsum(qbig));
% %     xlabel('r [Mpc]')    
% %     ylabel('q [M_\odot*(km/sec)^2]')
% %     title(sprintf('%s Cumulative Cooling Profile',clname)); 
% %     saveas(gcf,sprintf('%s/%s_cumcoolprof.png',outpath,clname));
% %     
%     coolprofs{i,1}=clname;
%     coolprofs{i,2}=squeeze(qnr(1,:,:));
%     coolprofs{i,3}=squeeze(qnr(2,:,:));
%     coolprofs{i,4}=squeeze(qnr(3,:,:));
%     coolprofs{i,5}=squeeze(qnr(4,:,:));
%     coolprofs{i,6}=qnrb;
%     
%     clear q1 q2 q4 q8 qbig rq
%     close gcf
%     
%     
% end
% save(sprintf('mat_files/cooling_profiles_%s.mat',aexp),'coolprofs');
