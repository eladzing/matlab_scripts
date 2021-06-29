%function  test_prof(cluster)



cluster='CL101';
type='csf';
aexp='a1';

new_env(cluster,type,aexp);

rp=0:0.0001:10;

[rog rot] =read_RHO_Profiles(rp);
tp = read_T_Profile(rp);

rott=rot(rot>0);
rpp=rp(rot>0);


%tpp=tp(tp>0);
figure;
loglog(rpp,rott)
%figure;
%loglog(rp,tp)

