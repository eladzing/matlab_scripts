dd=20;
rv=get_rvir

for i=1:dd
r=2*rv*i/dd
[VcmX VcmY VcmZ]=Vcm_full([0 r],[0 0 0],'nohub');
vv(end+1,:)=[VcmX VcmY VcmZ];
end