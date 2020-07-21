list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

%list2=[103 105 106 107];

result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout';
%result_dir='/home/eladzing/work/sshfs/titan3/cold_flows/tit3/printout/';
%result_dir='/home/titan3/eladzing/cold_flows/printout';

plotflag='noplot';
pflag='noprint';
%hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    

for id=1:length(list1)
    %halopath=sprintf(hpath,list1(id));
    clustername=sprintf('CL%d',list1(id))
    new_env(clustername);    
    stack{id,1}=clustername;
    switch list1(id)
        case {101,102,103,104,105,106,107,5}
            smallbox=2;bigbox=8;
        otherwise
            smallbox=1;bigbox=8;
    end
    time_profs
    stack{id,2}=RP/r200; %r
    stack{id,3}=tdyn1; %dynamical time in Gyr
    stack{id,4}=tcool; %cooling time in Gyr
    stack{id,5}=tcon; %contraction time in Gyr
    stack{id,6}=tdyn3; %r/cs  time in Gyr
    %%stack{id,5}=Tp./t200; %Temperature
    stack{id,7}=mdfull/Mg200; %full mass flux
    stack{id,8}=Mdp/Mg200; %approximate raidative mass flux
    %close all
end
save('mat_files/stk_times.mat','stack');

  
