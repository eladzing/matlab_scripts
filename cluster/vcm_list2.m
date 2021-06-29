%list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
list=11;
%vf=[1.39 1.681 1.239 1.527 2.172 2.009 2.624 2.87 3.037 2.363 1.95 2.954];



nb=50;xl=0.626*4/0.7;
vr=xl/nb:xl/nb:xl;

% vcmlist array indices:
% 1: cluster id 
% 2: r/rv vx vy vz 
vcmlist2=zeros(length(list),4,length(vr));

for i=1:length(list);
    cl=sprintf('CL%d',list(i));
    
    new_env(cl,'csf','a1','win');
    
    rv=get_rvir;
    vcmlist2(i,1,:)=vr./rv;
    for j=1:length(vr)
        [VcmX VcmY VcmZ]=Vcm_full(vr(j),[0 0 0],'hub');
        vcmlist2(i,2:4,j)=[VcmX VcmY VcmZ];
        clear VcmX VcmY VcmZ
    end
    
    
end