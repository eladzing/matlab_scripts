list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout/';
set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);

%pflag='print';
hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    

%rmfac=0.0;rxfac=0.02;
%tmin=-0.6;tmax=0.6;
%len=1000;
%tmin=-2;tmax=0.6;

stack={};

for id=1:length(list1)
    halopath=sprintf(hpath,list1(id));
    clustername=sprintf('CL%d',list1(id))
    stack{id,1}=clustername;
    switch list1(id)
        case {101,102,103,104,105,106,107,5}
            smallbox=2;
        otherwise
            smallbox=1;
    end
    [flprof rprof m200 mg200 r200]= flux_profile_allbox(halopath,smallbox);
    stack{id,2}=rprof./r200;
    stack{id,3}=flprof./mg200.*1e9; %%units are in 1/Gyr
end

%% plot the lines together
%figure;

%for i=1:size(stack,1)
%    if i==2
%        hold on
%    end
%    semilogx(stack{i,2},stack{i,3});
%end   
%hold off
  
%% find mean

%lbin=-3:0.05:log10(3);
%RP=10.^lbin;
%stx=zeros(size(RP));
%for i=1:N
%    stx=stx+spline(stack{i,2},stack{i,3},RP);
%end
%stx=stx./N;
