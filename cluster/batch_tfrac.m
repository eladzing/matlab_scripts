list1=[101 102 103 104 105 106 107 3 6 7 9 10 11 14 24];

%list2=[103 105 106 107];

result_dir='/home/alf/eladzing/tmp_results';

plotflag='plot';
pflag='print';
hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    

rmfac=0.0;rxfac=2;
tmfac=-0.5;

for id=1:length(list1)
    halopath=sprintf(hpath,list1(id));
    clustername=sprintf('CL%d',list1(id))
    stack{id,1}=clustername;
    switch list1(id)
        case {101,102,103,104,105,106,107,5}
            smallbox=2;bigbox=8;
        otherwise
            smallbox=1;bigbox=8;
    end
    tfrac_prof
    stack{id,2}=rbin;
    stack{id,3}=tfprof2./mnorm2;
end


  