list=[101 103 105 106 3 5 6 9 10 11 14 24];

vf=[1.39 1.681 1.239 1.527 2.172 2.009 2.624 2.87 3.037 2.363 1.95 2.954];


%nb=50;xl=4/0.7;
%vr=xl/nb:xl/nb:xl;

% vcmlist array indices:
% 1: cluster id 
% 2: r/rv vx vy vz 
bl=8/0.7/2;
vcmlist3=zeros(length(list),5,4);

for i=1:length(list);
    cl=sprintf('CL%d',list(i));
    
    new_env(cl,'csf','a1');
    
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