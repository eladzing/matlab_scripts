%list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%flimx=[10 7 4.8 1 5.37 1.44 1.15 2 0.6 4.8 0.3 0.3 0.03 2 0.23 0.06];
%flimn=[25 14.4 10 2.5 2.88 4.8 1.15 1.15 0.5 2.5 0.3 0.73 0.06 9.6 0.15 0.06];
%list2=[103 105 106 107];
%global CLUSTER;
%CLUSTER='CL%s';
result_dir='/home/titan3/eladzing/cold_flows/printout';
stack_dir='/home/carina/eladzing_a/cold_flows/datacube/matlab/mat_files';

aexp='a1';typ='csf';
hubbleflag='hub';

pflag='noprint';
hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    

%rmf=[0 0.02 0.1 0.6];
%rxf=[0.01 0.1 0.3 1.2];

rmf=[0 0.02 0.35];
rxf=[0.01 0.25 1.2];

ptags={'core_clean_vzone2','clean_vzone2','outer_clean_vzone2'};

%% zone 1 : core 
%%rmfac=0;rxfac=0.01;
%%printag='core_clean';
%% zone 2 : inner 
%rmfac=0.02;rxfac=0.1;
%printag='inner_clean';
%% zone 2a : inner 
%rmfac=0.1;rxfac=0.3;
%printag='inner_clean2';
%% zone 3 : outer 
%rmfac=0.6;rxfac=1.2;
%printag='outer_clean';



%fln=-0.16e-13;flx=0.16e-13;
fln=-0.08;flx=0.08;
tn=-1.7;tx=0.8;
Rvcm=[0 0.2]; %Radii fractions for Vcm calculation 

for ii=1:3
    rmfac=rmf(ii);
    rxfac=rxf(ii);
    printag=char(ptags(ii));  
        
%%for id=1:length(list1)
    %halopath=sprintf(hpath,list1(id));
    CLST='CL6';
    %%CLUSTER=sprintf(CLUSTER,num2str(list1(id)))
    new_env(CLST,typ,aexp);
    
    %stack{id,1}=CLUSTER;
    %switch list1(id)
    %    case {101,102,103,104,105,106,107,5}
    %        smallbox=2;bigbox=8;
    %    otherwise
            smallbox=1;bigbox=8;
    %end
    %fln=-1.*flimn(id);flx=flimx(id);
    flux_birds_2   
    %stack{id,2}=tax;
    %stack{id,3}=mdhist;
    %close all
    
    
end
%save(sprintf('%s/stk_mdhist_%s.mat',stack_dir,printag),'stack');
%end
  