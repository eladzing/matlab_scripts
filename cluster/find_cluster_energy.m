%% find cluster energy
% find the energy profile of a cluster. 


list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
mask =[ 0   0   0   0   0   0   0 0 0 1 0 0  0  0  0  0];
ll=sum(mask);
energyprofs=cell(ll,11);
   
%% global defineitions 
typ='csf';
aexp='a1';

hubflag='hub';
vcenterflag='defualt';
center=[0,0,0];


kb=1.38066e-16; %erg/K
mu=0.5926; %mass per particle for primordial composition
mp=1.672649e-24; % gram
yr=3.155815e7; %  sec 
pc=3.0856e18; % center
km=1e5; %center
Ms=1.989e33; %gr
xi=0.5185185; 

%%numerical factors converting to units of M_sun*(km/sec)^2 per Mpc^3
f_th=1.5.*kb./(mu.*mp)/km.^2 ;%E_th factor converts to units cited above
f_pv=kb./(mu.*mp.*km.*1.0e6.*pc).*yr; %units of M_sun*(km/sec)^2 per Mpc^3 per yr
fu=km./(1e6.*pc).*yr;
f_cl=(xi./(mu.*mp)).^2*yr*Ms/(1e6*pc)^3/km^2*1.e-22;

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

    for i=1:4
    
        boxx=2.^(i-1)
        
        vol=(boxx./NCELL./hub).^3; %in Mpc^3 proper
        
        rp(i,:)=(ri-0.5).*(boxx./2./NCELL./hub);

        %% find potential energies 
        [pdm pgs pdot ro]=poteng(boxx);
        epg=(pgs+pdm).*ro.*vol;
        edot=ro.*pdot.*vol;
        
        epot(i,:)=MAKE_PROFILE_FROM_CUBE(epg);
        epdot(i,:)=MAKE_PROFILE_FROM_CUBE(edot);

        clear epg edot pgs pdmpdot
        %% find kinetic energy 
        
        % calculate the spherical velocity components
        [Vr Vtheta Vphi] = V_decomposition(boxx,vcenterflag,hubflag,center);
        
        % transfer to spherical coordinates 
%         Vr=cart2sphere(Vr);
%         Vtheta=cart2sphere(Vtheta);
%         Vphi=cart2sphere(Vphi);
%         
%         calculate the velocity dispersion        
%         [sigt rmst vtb]=find_dispersion(vt(1:end-1,:,:),ros);
%         

        ekr=0.5.*ro.*vol.*(Vr.^2+Vtheta.^2+Vphi.^2);
        ekt=0.5.*ro.*vol.*(Vtheta.^2+Vphi.^2);
        
        % make profile
        kpr(i,:)=MAKE_PROFILE_FROM_CUBE(ekr);
        kpt(i,:)=MAKE_PROFILE_FROM_CUBE(ekt);
        
        clear ekr ekt Vr Vtheta Vphi
        
        %% calculate thermal 
        Tp=T(boxx);
        eth=ro.*Tp.*vol.*f_th;
        thp(i,:)=MAKE_PROFILE_FROM_CUBE(eth);
        
        % Pdv part

        Pdv=-1.0.*ro.*Tp.*divV(boxx).*f_pv.*vol;
        ppr(i,:)=MAKE_PROFILE_FROM_CUBE(Pdv);
        
        
        %% calculate mass profile
        mprof(i,:)=MAKE_PROFILE_FROM_CUBE(ro.*vol);
        clear ro Tp vol eth Pdv
    end
    
    load(sprintf('%s/virial%d_%s', HALO_PATH,8,aexp));
    [rq thb]=concat_qprofs(thp(1,:),thp(2,:),thp(3,:),thp(4,:));
    [rr pvrate]=concat_qprofs(ppr(1,:),ppr(2,:),ppr(3,:),ppr(4,:));
    
    [rr ekrbig]=concat_qprofs(kpr(1,:),kpr(2,:),kpr(3,:),kpr(4,:));
    [rr ektbig]=concat_qprofs(kpt(1,:),kpt(2,:),kpt(3,:),kpt(4,:));
    [rr mbig]=concat_qprofs(mprof(1,:),mprof(2,:),mprof(3,:),mprof(4,:));

    ekprof(:,1)=ekrbig;
    ekprof(:,2)=ektbig;
    
    
    
    [rr epbig]=concat_qprofs(epot(1,:),epot(2,:),epot(3,:),epot(4,:));
    [rr edbig]=concat_qprofs(epdot(1,:),epdot(2,:),epdot(3,:),epdot(4,:));  
    
    %get cooling profile 
    load(sprintf('mat_files/cooling_profiles_%s.mat',aexp));
    qprof=coolprofs{id,6}; % 1st index rprofile, second indexd is cooling;
    
    explain='contains: 2-name;3-rporifile;4-kinetic(r,t);5-potential;6-thermal;7-cooling;8-pdv;9-epdot;10-virials'
    energyprofs{id,1}=explain;
    energyprofs{id,2}=cluster;
    energyprofs{id,3}=rq;
    energyprofs{id,4}=ekprof;
    energyprofs{id,5}=epbig;
    energyprofs{id,6}=thb;
    energyprofs{id,7}=qprof(2,:);
    energyprofs{id,8}=pvrate;
    energyprofs{id,9}=edbig;
    energyprofs{id,10}=mbig;
    energyprofs{id,11}=[MVIR,RVIR,VVIR,TVIR];
        
    
end