%list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

%list1=[101 3  6  11 14];

list1=[101];% 3  6  11 14];

load(sprintf('mat_files/cooling_profiles_%s.mat','a1'));

typ='csf';
aexp='a1';

hubflag='hub';
vcmflag='vcm';

smallbox=1;
bigbox=8;
cen_fac=0.05;

%%some global definitions
kb=1.38066e-16; %erg/K
mu=0.5926; %mass per particle for primordial composition
mp=1.672649e-24; % gram
yr=3.155815e7; %  sec 
pc=3.0856e18; % cm
km=1e5; %cm
Ms=1.989e33; %gr
xi=0.5185185; 

%%numerical factors converting to units of M_sun*(km/sec)^2 per yr
f_th=1.5.*kb./(mu.*mp); %thermal energy factor
fu=km./(1e6.*pc).*yr;
f_cl=(xi./(mu.*mp)).^2*yr*Ms/(1e6*pc)^3/km^2*1.e-22;

cm=[0,0,0];

for id=1:length(list1)
    cluster=sprintf('CL%d',list1(id))   
    eng_stack{id,1}=cluster;
    new_env(cluster,typ,aexp);

    qb=coolprofs{id,6};
    rq=qb(1,:);
    q=qb(2,:);
    clear qb;
    [inedg rinedg]=find_inner_cooling_edge(rq,q);
    
    global HALO_PATH
    global aexpn;
    global hub;
    global NCELL;   
    
    [qb1 qb2 qb4 qb8] = cooling_profs(cluster,typ,aexp,cm,cen_fac) ;
    %if ~exist('rvcm','var')
    %    rvcm=[0 1];
    %end

    if ~exist('hubflag','var')
        hf=0;
    else
    hf=(strcmpi(hubflag,'hub'));    
    end

    Ekr=[];
    Ekt=[];
    Ept=[];
    Eth=[];
    qcl=[];
    qcl1=[];
    
    boxx=smallbox

    if strcmp(vcmflag,'vcm') %%find Vrcm within Rvir
        %rvc=rvcm.*get_rvir();
        vrcm=Vrcm_rvir_SPHERE(boxx);
    else
        vrcm=zeros(size(full_ro));
    end
    %clear rvc;

    load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));

    % load relevant cubes
    ds=ds_sphere(boxx);
    vr=Vr_sphere(boxx)+hf.*hubble_sphere(boxx)-vrcm;
    vrds=vr.*ds;
    drds=ds.*boxx/hub/2/NCELL;
    clear ds vr 

    rosp=RHOG_sphere(boxx);
    Tsp=T_sphere(boxx);

    [Eks Ekrs] = Ekin_sphere(boxx);
    %Eks=Eks %.*rosp.*vrds;
    %Ekrs=Ekrs%.*rosp.*vrds;
    Ekts=(Eks-Ekrs).*vrds; clear Eks
    Ekrs=Ekrs.*vrds;

    eths=rosp.*Tsp.*vrds;
    [rpr po] = potential_prof(R_Profile);

    flx=rosp.*vrds; %%flux 
    
    zmet=ZIa_sphere(boxx)+ZII_sphere(boxx);
    tt=log10(Tsp);
    rosq=rosp.*rosp.*drds;    
    
    clear vrds Tsp rosp drds

    for ridx = 1:(length(R_Profile)-1)
        %% kinetic terms
        Ekr(end+1)=sum(sum(squeeze(Ekrs(ridx,:,:))));
        Ekt(end+1)=sum(sum(squeeze(Ekts(ridx,:,:))));
    
        %% Thermal term
        Eth(end+1)=sum(sum(squeeze(eths(ridx,:,:))));

        %% Potential terms
        pot=interp1(rpr,po,R_Profile(ridx),'spline');
        Ept(end+1)=sum(sum(squeeze(flx(ridx,:,:).*pot)));
    
        %% Cooling
        if R_Profile(ridx)>=cen_fac*get_rvir  
            rom=squeeze(rosq(ridx,:,:));
            zm=squeeze(zmet(ridx,:,:));
            tm=squeeze(tt(ridx,:,:));
            lamb=interp2(zlam,tlam,lambda,zm,tm,'spline');
            mask=(tm>4.0);
            lamb2=10.^(lamb+22).*mask;
            qcl1(end+1)=sum(sum(lamb2.*rom));
        else
            qcl1(end+1)=0;
        end
            
    
    end
    rp=R_Profile(1:(length(R_Profile)-1));
    rp1=rp;
    qcl=qb1(1:(length(qb1)-1));
    
    for ii=(log2(smallbox)+1):log2(bigbox)  
        boxx=2^ii
    
        load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));
    
        % load relevant cubes
        ds=ds_sphere(boxx);
        vr=Vr_sphere(boxx)+hf.*hubble_sphere(boxx)-vrcm;
        vrds=vr.*ds;
        drds=ds.*boxx/hub/2/NCELL;
        clear ds vr

        rosp=RHOG_sphere(boxx);
        Tsp=T_sphere(boxx);

        [Eks Ekrs] = Ekin_sphere(boxx);
        Ekts=(Eks-Ekrs).*rosp.*vrds; clear Eks
        Ekrs=Ekrs.*rosp.*vrds;
    
        eths=rosp.*Tsp.*vrds;
        [rpr po] = potential_prof(R_Profile);

        flx=rosp.*vrds; 

        zmet=ZIa_sphere(boxx)+ZII_sphere(boxx);
        tt=log10(Tsp);
        rosq=rosp.*rosp.*drds;    
        
        clear rosp vrds Tsp drds

    
        for ridx =(length(R_Profile)/2):(length(R_Profile)-1)
            %% kinetic terms
            Ekr(end+1)=sum(sum(squeeze(Ekrs(ridx,:,:))));
            Ekt(end+1)=sum(sum(squeeze(Ekts(ridx,:,:))));

            %% Thermal term
            Eth(end+1)=sum(sum(squeeze(eths(ridx,:,:))));

            %% Potential terms
            pot=interp1(rpr,po,R_Profile(ridx),'spline');
            Ept(end+1)=sum(sum(squeeze(flx(ridx,:,:).*pot)));
    
            %% Cooling
            if R_Profile(ridx)>=cen_fac*get_rvir  
                rom=squeeze(rosq(ridx,:,:));
                zm=squeeze(zmet(ridx,:,:));
                tm=squeeze(tt(ridx,:,:));
                lamb=interp2(zlam,tlam,lambda,zm,tm,'spline');
                mask=(tm>4.0);
                lamb2=10.^(lamb+22).*mask;
                qcl1(end+1)=sum(sum(lamb2.*rom));
            else
                qcl1(end+1)=0;
            end
        end
        rp=cat(2,rp,R_Profile(length(R_Profile)/2:length(R_Profile)-1));
        switch ii    
            case 1
                rp2=R_Profile;
                qb=qb2;
            case 2
                rp4=R_Profile;
                qb=qb4;
            case 3
                rp8=R_Profile;
                qb=qb8;
        end
        qcl=cat(2,qcl,qb(length(qb)/2:length(qb)-1));
    end

    if bigbox<8
        load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexpn));
    end

    Ekr=Ekr.*fu;
    Ekt=Ekt.*fu;
    Ept=Ept.*fu;
    Eth=Eth.*f_th.*fu./km.^2;
    qcl1=qcl1.*f_cl;
    
    eng_stack{id,2}=rp/RVIR;
    eng_stack{id,3}=Ekr;
    eng_stack{id,4}=Ekt;
    eng_stack{id,5}=Ept;
    eng_stack{id,6}=Eth;
    eng_stack{id,7}=qcl;
    eng_stack{id,8}=[MVIR,RVIR,VVIR,TVIR];
    
%    clear rp Ekr Ekt Ept Eth qcl qb qb1 qb2 qb4 qb8
end
%save(sprintf('mat_files/stk_eng_profs_%s.mat',aexp),'eng_stack');