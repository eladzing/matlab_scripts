list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

vf=[1.5 1.8 1.3 1.2 1.08 1.5 1.5 1.3 1.5 1.24 1.39 1.7 1.8 1.5 2.3 1.5];
    


%nb=50;xl=4/0.7;
%vr=xl/nb:xl/nb:xl;

% vcmlist array indices:
% 1: cluster id 
% 2: r/rv vx vy vz 
bl=0.626*8/0.7/2;
vcmlist3=zeros(length(list),5,4);

for i=1:length(list);
    cl=sprintf('CL%d',list(i));
    
    new_env(cl,'csf','a06');
    global zred
    
    rv=get_rvir;
    vr=[vf(i) 1 2 3].*rv;
    vr(vr>bl)=bl;
    vcmlist3(i,1,:)=vr;
    for j=1:4
        
        [VcmX VcmY VcmZ]=Vcm_full(vr(j),[0 0 0],'hub');
        vcmlist3(i,2:4,j)=[VcmX VcmY VcmZ];
        vcmlist3(i,5,j)=sqrt(VcmX.^2+VcmY.^2+VcmZ.^2);
        clear VcmX VcmY VcmZ
    end
    
    
end