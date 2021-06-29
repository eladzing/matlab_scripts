%% find cluster energy
% find the energy budget of a eUlarian shell


list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
mask =[ 0   0   0   0   0   0   0  0 0 1 0 0  0  0  0  0];
ll=length(list); 
energyprofs=cell(ll,10);
   
%% global defineitions 
typ='csf';
aexp='a1';

hubflag='hub';
vcenterflag='defualt';
center=[0,0,0];

units;

% in units of Rv 
rsh=[0 0.5];


for id=1:length(list)
    if(~mask(id)) 
        continue
    end
        
    cluster=sprintf('CL%d',list(id))   

    %begin calculation 
    new_env(cluster,typ,aexp);

    global HALO_PATH
    %global aexpn;
    global hub;
    global NCELL;   
    ri=1:NCELL;
          
    % find box indices for shells
    rshell=rsh.*get_rvir;
    ib=8*ones(size(rshell));
    ib(rshell<2./hub)=4;
    ib(rshell<1./hub)=2;
    ib(rshell<0.5./hub)=1;
    
    %check inner - 
    rri=(ri-0.5).*(ib(1)./2./NCELL./hub);
    r_ind(1)=1; %find(rri<rshell(1),1,'last');
    rro=(ri-0.5).*(ib(2)./2./NCELL./hub);
    r_ind(2)=find(rro>rshell(2),1,'first');
    
    r_ind(r_ind==NCELL)=NCELL-1;
    rshell=[rri(r_ind(1)) rro(r_ind(2))]       
    
    
    %ib=ceil(2.^ceil(log2(rshell/0.5))); 
    %ib(ib==0)=1; % fix for zero value
    
    %get cooling profile and find total cooling in shell  
    load(sprintf('mat_files/cooling_profiles_%s.mat',aexp));
    qprof=coolprofs{id,6}; % 1st index rprofile, second indexd is cooling;
    rpn=qprof(1,:)./hub;qp=qprof(2,:);
    ind=find(rpn>=rshell(1) & rpn<=rshell(2));
    qtot=sum(qp(ind));

   
    
    for box=1:4;
        boxx=2^(box-1);
        disp(boxx);
        [pdm pgs pdot ro]=poteng(boxx);
        clear pdot             
        
        % calculate the spherical velocity components
        [Vr Vtheta Vphi] = V_decomposition(boxx,vcenterflag,hubflag,center);
        
        vsq=Vr.^2+Vtheta.^2 + Vphi.^2;
        
        clear Vtheta Vphi
        
        tp=T(boxx);
         
        Pr=ro.*tp.*f_pr;
        eint=ro.*tp.*f_th;
        
        clear tp 

        ecube=eint+0.5.*ro.*vsq+ro.*(pdm+pgs)+Pr;
        
       clear Pr eint ro vsq pdm pgs 
         
        esphere=cart2sphere(ecube);
        vrs=cart2sphere(Vr);
        
        %eprof(box,:)=MAKE_PROFILE_FROM_CUBE(esphere.*vrs.*ds_sphere(boxx));
        
        eprof(box,:)=sum(sum(esphere.*vrs.*ds_sphere(boxx),2),3);
        
        switch(boxx)
            case(ib(1))
                einner=eprof(r_ind(1));
            case(ib(2))
                eouter=eprof(r_ind(2));
        end
        
        clear esphere vrs ecube        
        
         
    end 
    [rr epbig]=concat_qprofs(eprof(1,:),eprof(2,:),eprof(3,:),eprof(4,:));
     
         % transfer to spherical coordinates 
         %Vr=cart2sphere(Vr);
         
 
         
        
       
%         
%         get cooling 
%         
%         
%         
%     end
%     
%     
%     ri=1:NCELL;
% 
%     for i=1:4
%     
%         boxx=2.^(i-1)
%         
%         vol=(boxx./NCELL./hub).^3; %in Mpc^3 proper
%         
%         rp(i,:)=(ri-0.5).*(boxx./2./NCELL./hub);
% 
%         %% find potential energies 
%         [pdm pgs pdot ro]=poteng(boxx);
%         epg=(pgs+pdm).*ro.*vol;
%         edot=ro.*pdot.*vol;
%         
%         epot(i,:)=MAKE_PROFILE_FROM_CUBE(epg);
%         epdot(i,:)=MAKE_PROFILE_FROM_CUBE(edot);
% 
%         clear epg edot pgs pdmpdot
%         %% find kinetic energy 
%         
%         % calculate the spherical velocity components
%         [Vr Vtheta Vphi] = V_decomposition(boxx,vcenterflag,hubflag,center);
%         
%         % transfer to spherical coordinates 
% %         Vr=cart2sphere(Vr);
% %         Vtheta=cart2sphere(Vtheta);
% %         Vphi=cart2sphere(Vphi);
% %         
% %         calculate the velocity dispersion        
% %         [sigt rmst vtb]=find_dispersion(vt(1:end-1,:,:),ros);
% %         
% 
%         ekr=0.5.*ro.*vol.*(Vr.^2+Vtheta.^2+Vphi.^2);
%         ekt=0.5.*ro.*vol.*(Vtheta.^2+Vphi.^2);
%         
%         % make profile
%         kpr(i,:)=MAKE_PROFILE_FROM_CUBE(ekr);
%         kpt(i,:)=MAKE_PROFILE_FROM_CUBE(ekt);
%         
%         clear ekr ekt Vr Vtheta Vphi
%         
%         %% calculate thermal 
%         Tp=T(boxx);
%         eth=ro.*Tp.*vol.*f_th;
%         thp(i,:)=MAKE_PROFILE_FROM_CUBE(eth);
%         
%         % Pdv part
% 
%         Pdv=-1.0.*ro.*Tp.*divV(boxx).*f_pv.*vol;
%         ppr(i,:)=MAKE_PROFILE_FROM_CUBE(Pdv);
%         
%         clear ro Tp vol eth Pdv
%     end
%     
%     load(sprintf('%s/virial%d_%s', HALO_PATH,8,aexp));
%     [rq thb]=concat_qprofs(thp(1,:),thp(2,:),thp(3,:),thp(4,:));
%     [rr pvrate]=concat_qprofs(ppr(1,:),ppr(2,:),ppr(3,:),ppr(4,:));
%     
%     [rr ekrbig]=concat_qprofs(kpr(1,:),kpr(2,:),kpr(3,:),kpr(4,:));
%     [rr ektbig]=concat_qprofs(kpt(1,:),kpt(2,:),kpt(3,:),kpt(4,:));
% 
%     ekprof(:,1)=ekrbig;
%     ekprof(:,2)=ekrbig;
%     
%     
%     
%     [rr epbig]=concat_qprofs(epot(1,:),epot(2,:),epot(3,:),epot(4,:));
%     [rr edbig]=concat_qprofs(epdot(1,:),epdot(2,:),epdot(3,:),epdot(4,:));  
%     
%     %get cooling profile 
%     load(sprintf('mat_files/cooling_profiles_%s.mat',aexp));
%     qprof=coolprofs{id6}; % 1st index rprofile, second indexd is cooling;
%     
%     explain='contains: 2-name;3-rporifile;4-kinetic(r,t);5-potential;6-thermal;7-cooling;8-pdv;9-epdot;10-virials'
%     energyprofs{id,1}=explain;
%     energyprofs{id,2}=cluster;
%     energyprofs{id,3}=rq;
%     energyprofs{id,4}=ekprof;
%     energyprofs{id,5}=epbig;
%     energyprofs{id,6}=thb;
%     energyprofs{id,7}=qprof(2,:);
%     energyprofs{id,8}=pvrate;
%     energyprofs{id,9}=edbig;
%     energyprofs{id,10}=[MVIR,RVIR,VVIR,TVIR];
        
    
 end