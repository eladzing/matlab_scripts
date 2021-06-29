rc=mk_rcube(8,ones([256 256 256]));
per=[];
nn1=[];nn2=[];
acs=[];
for ii=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

    clustername=sprintf('CL%d',ii)
    new_env(clustername);  
    
    rv=get_rvir();
    
    rv_ind=find(rc<=rv);
    rv_ind2=find(rc<=0.2.*rv);
    
    [VcmX VcmY VcmZ] = V_Vcm_r(8,rv_ind); 
    [VcmX2 VcmY2 VcmZ2] = V_Vcm_r(8,rv_ind2);
    
    n1=norm([VcmX VcmY VcmZ]);
    n2=norm([VcmX2 VcmY2 VcmZ2]);
    nn1(end+1)=n1;nn2(end+1)=n2;
    per(end+1)=100*abs(1-n1./n2);
    
    acs(end+1)=acos([VcmX VcmY VcmZ]*[VcmX2 VcmY2 VcmZ2]'./(n1*n2))*180./pi;
    
end