

halopath='/home/alf/eladzing/data/kravtsov/clusters/%s';

[flprof roprof vrprof r_prof rvir mvir]= full_profs(sprintf(halopath,cluster),smallbox);


rvir

vc=[0.02 0.1].*rvir

ror=spline(r_prof,roprof,vc)
vrr=spline(r_prof,vrprof,vc)
flr=spline(r_prof,flprof,vc)


MGAS=read_MGAS_Profile(sprintf(halopath,cluster),rvir)
mvir
m_norm=-560*(mvir./1e13).^(0.15).*(MGAS./1e13)  %virial accretion rate in Msun/yr



