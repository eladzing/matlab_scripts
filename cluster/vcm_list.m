list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
vr=0.1:0.1:3.0;

vcmlist=zeros(length(list),3,length(vr));

for i=1:length(list);
    cl=sprintf('CL%d',list(i));
    
    new_env(cl,'csf','a06');
    
    rv=get_rvir;
    
    for j=1:length(vr)
        [VcmX VcmY VcmZ]=Vcm_full(vr(j).*rv,[0 0 0],'hub');
        vcmlist(i,:,j)=[VcmX VcmY VcmZ];
        clear VcmX VcmY VcmZ
    end
    
    
end