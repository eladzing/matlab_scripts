list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

%list2=[103 105 106 107];

result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout/';
%result_dir='/home/eladzing/work/sshfs/titan3/cold_flows/tit3/printout/';

plotflag='plot';
pflag='print';
%hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    

%rmfac=0.0;rxfac=2;
%tmfac=-0.5;

rmf=[0 0.02 0.1 0.6];
rxf=[0.01 0.1 0.3 1.2];

ptags={'core_clean','inner_clean','inner_clean2','outer_clean'};

%% zone 1 : core 
%rmfac=0;rxfac=0.01;
%printag='core_clean';
%% zone 2 : inner 
%rmfac=0.02;rxfac=0.1;
%printag='inner_clean';
%% zone 2a : inner 
%rmfac=0.1;rxfac=0.3;
%printag='inner_clean2';
%% zone 3 : outer 
%rmfac=0.6;rxfac=1.2;
%printag='outer_clean';

for ii=1:4
    
    rmfac=rmf(ii)
    rxfac=rxf(ii)
    printag=char(ptags(ii))    

for id=1:length(list1)
    %halopath=sprintf(hpath,list1(id));
    clustername=sprintf('CL%d',list1(id))
    new_env(clustername);    
    %stack{id,1}=clustername;
    switch list1(id)
        case {101,102,103,104,105,106,107,5}
            smallbox=2;bigbox=8;
        otherwise
            smallbox=1;bigbox=8;
    end
    time_bird
    %stack{id,2}=RP/r200; %r
    %stack{id,3}=tdyn1; %dynamical time in Gyr
    %stack{id,4}=tcool; %cooling time in Gyr
    %stack{id,5}=tcon; %contraction time in Gyr
    %stack{id,5}=Tp./t200; %Temperature
    %stack{id,6}=mdfull/Mg200; %full mass flux
    %stack{id,7}=Mdp/Mg200; %approximate raidative mass flux
    close all
end
end

  
