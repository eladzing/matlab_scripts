%list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

%ist1=[101 3  6  11 14];

list1=[6];% 3  6  11 14];

%load(sprintf('mat_files/cooling_profiles_%s.mat','a1'));

typ='csf';
aexp='a1';

hubflag='hub';
vcmflag='vcm';

%cen_fac=0.05;%%

%%some global definitions
kb=1.38066e-16; %erg/K
mu=0.5926; %mass per particle for primordial composition
mp=1.672649e-24; % gram
yr=3.155815e7; %  sec 
pc=3.0856e18; % cm
km=1e5; %cm
Ms=1.989e33; %gr
xi=0.5185185; 

%%numerical factors converting to units of M_sun*(km/sec)^2 per Mpc^3
%f_th=1.5.*kb.*Ms./(mu.*mp)    ; %E_th factor converts to erg per Mpc^3
f_th=1.5.*kb./(mu.*mp)/km.^2 ;%E_th factor converts to units cited above
f_pv=kb./(mu.*mp.*km.*1.0e6.*pc).*yr; %units of M_sun*(km/sec)^2 per Mpc^3 per yr
fu=km./(1e6.*pc).*yr;
f_cl=(xi./(mu.*mp)).^2*yr*Ms/(1e6*pc)^3/km^2*1.e-22;

cm=[0,0,0];

for id=1:length(list1)
    cluster=sprintf('CL%d',list1(id))   
    
    potential_eng_stack{id,1}=cluster;
    
    new_env(cluster,typ,aexp);
   
    
    global HALO_PATH
    global aexpn;
    global hub;
    global NCELL;   
        

    
    ri=1:NCELL;
    epot=[];rp=[];edot=[];
    for i=1:4
            
        boxx=2.^(i-1)
        
        vol=(boxx./NCELL./hub).^3; %in Mpc^3 proper
        
        rp(i,:)=(ri-0.5).*(boxx./2./NCELL./hub);
        
        %find potential energies 
        [pdm pgs pdot ro]=poteng(boxx);
        epg=(pgs+pdm).*ro.*vol;
        edot=ro.*pdot.*vol;
    
        epot(i,:)=MAKE_PROFILE_FROM_CUBE(epg);
        epdot(i,:)=MAKE_PROFILE_FROM_CUBE(edot);
        
        %% kinetic part
        %[ek ekr]=Ekin(boxx);  %.*vol;
        %ek=ek.*vol;
        %ekr=ekr.*vol;
        %ekt=ek-ekr;
        %kpr(i,:)=MAKE_PROFILE_FROM_CUBE(ekr);
        %kpt(i,:)=MAKE_PROFILE_FROM_CUBE(ekt);
%        %kp(i,:)=MAKE_PROFILE_FROM_CUBE(ek);
        %clear ek ekt ekr
        
%         %% thermal part
%         ro=RHOG(boxx);
%         Tp=T(boxx);
%         
%         eth=ro.*Tp.*vol.*f_th;
%         thp(i,:)=MAKE_PROFILE_FROM_CUBE(eth);
%                     
%         %% Pdv part
%          %divu=divV(boxx);
%          Pdv=-1.0.*ro.*Tp.*divV(boxx).*f_pv.*vol;
%          ppr(i,:)=MAKE_PROFILE_FROM_CUBE(Pdv);
        
        clear ro epg edot vol pgs pdm pdot
    
    end
    
    [rq epbig]=concat_qprofs(epot(1,:),epot(2,:),epot(3,:),epot(4,:));
    [rq1 edbig]=concat_qprofs(epdot(1,:),epdot(2,:),epdot(3,:),epdot(4,:));
        
    rq=rq./hub;
    
    load(sprintf('%s/virial%d_%s', HALO_PATH,8,aexpn));
     
          
     potential_eng_stack{id,2}=rp;  %radii (in Mpc proper)
     potential_eng_stack{id,3}=epot; %potential energy
     potential_eng_stack{id,4}=epdot; % potential energy change rate
     potential_eng_stack{id,5}=rq; % big radii array
     potential_eng_stack{id,6}=epbig; %potential big profile
     potential_eng_stack{id,7}=edbig; %potential big profile
     potential_eng_stack{id,8}=[MVIR,RVIR,VVIR,TVIR];
     clear rp epbig edbig rq rp epot edot
       
end

%save(sprint.*f('mat_files/stk_eng_test_%s.mat',aexp),'eng_stack');