outpath='~/work/clusters/printout';

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
ll=length(list);
coolprofs=cell(ll,7);
typ='csf';
aexp='a1';
cm=[0,0,0];

nc=256;
ri=1:nc;
r1=(ri-0.5).*(1./2./nc);
r2=(ri-0.5).*(2./2./nc);
r4=(ri-0.5).*(4./2./nc);
r8=(ri-0.5).*(8./2./nc);


for i=1:length(list)
    qnr=zeros(4,2,nc);
    clname=sprintf('CL%d',list(i))
    [q1 q2 q4 q8] = cooling_profs(clname,typ,aexp,cm,0);   
   
    
    qnr(1,1,:)=r1;qnr(1,2,:)=q1;
    qnr(2,1,:)=r2;qnr(2,2,:)=q2;
    qnr(3,1,:)=r4;qnr(3,2,:)=q4;
    qnr(4,1,:)=r8;qnr(4,2,:)=q8;
    
    [rq qbig] = concat_qprofs(q1,q2,q4,q8);
    qnrb(1,:)=rq;qnrb(2,:)=qbig;
    [inedg rinedg]=find_inner_cooling_edge(rq,qbig,clname)
    
    figure;
    loglog(r1,q1,r2,q2,r4,q4,r8,q8)
    hold on
        loglog([rinedg rinedg],[min(q8) max(q8)],'--r');
    hold off
    xlabel('r [Mpc]')    
    ylabel('q [M_\odot*(km/sec)^2]')
    title(sprintf('%s Cooling Profiles',clname)); 
    saveas(gcf,sprintf('%s/%s_coolprof.png',outpath,clname));
    
    figure;
    loglog(r1,q1,r2,q2./2,r4,q4./4,r8,q8./8)
    hold on
        loglog([rinedg rinedg],[min(q8) max(q8)],'--r');
    hold off
    xlabel('r [Mpc]')    
    ylabel('q [M_\odot*(km/sec)^2]')
    title(sprintf('%s Cooling Profiles normalized',clname)); 
    saveas(gcf,sprintf('%s/%s_coolprof2.png',outpath,clname));

    
    
    %     figure;
%     loglog(rq,cumsum(qbig));
%     xlabel('r [Mpc]')    
%     ylabel('q [M_\odot*(km/sec)^2]')
%     title(sprintf('%s Cumulative Cooling Profile',clname)); 
%     saveas(gcf,sprintf('%s/%s_cumcoolprof.png',outpath,clname));
%     
    coolprofs{i,1}=clname;
    coolprofs{i,2}=squeeze(qnr(1,:,:));
    coolprofs{i,3}=squeeze(qnr(2,:,:));
    coolprofs{i,4}=squeeze(qnr(3,:,:));
    coolprofs{i,5}=squeeze(qnr(4,:,:));
    coolprofs{i,6}=qnrb;
    coolprofs{i,7}=[inedg rinedg];
    clear q1 q2 q4 q8 qbig rq inedg rinedg
    pause    
    close all
    
    
end
%save(sprintf('mat_files/cooling_profiles_%s.mat',aexp),'coolprofs');
