
load mat_files/stk_times.mat

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%      1   2   3   4   5   6   7  8 9 10111213 14 15 16   

cflist=[104 105 107 3 6 7 10 14;4 5 7 8 10 11 13 15];
ncflist=[101 102 103 106 5 9 11 24;1 2 3 6 9 12 14 16];
%ind2=RP/r200; %r
%ind3=tdyn1; %dynamical time in Gyr
%ind4=tcool; %cooling time in Gyr
%ind5=tcon; %contraction time in Gyr
%ind5=Tp./t200; %Temperature
%ind6=mdfull/Mg200; %full mass flux
%ind7=Mdp/Mg200; %approximate raidative mass flux


r_p=0.001:0.001:2.5;

for i=1:length(cflist);
    %clname=sprintf('CL%d',cflist(i));
    %for j=1:size(stack,1)
    %ind=find(strcmp(stack{:,1}
     
    tcool(i,:)=interp1(stack{i,2},stack{i,4},r_p,'spline');   
    tcon(i,:)=interp1(stack{i,2},stack{i,5},r_p,'spline');     
    tcool(i,:)=interp1(stack{i,2},stack{i,4},r_p,'spline');   
    tcon(i,:)=interp1(stack{i,2},stack{i,5},r_p,'spline'); 
