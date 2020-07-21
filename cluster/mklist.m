list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    

fac=[0.02 0.1];
%vc=zeros(length(list),2);
vc=[];

for id=1:length(list)  
    halopath=sprintf(hpath,list(id));
    load(sprintf('%s/virial%d', halopath, 8));
    vc(end+1,:)=fac.*(RVIR*1000);
end

