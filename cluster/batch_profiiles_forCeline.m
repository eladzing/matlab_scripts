rr=0.1:0.01:5;

list=[104,3,5,7,10,14];

for i=1:length(list)
    
    cl=list(i);
    new_env(cl);

[rog,~] = read_RHO_Profiles(rr);
tp=read_T_Profile(rr);

figure
loglog(rr.*1000,rog);
hold on 
loglog(1000.*get_rvir(500).*[1 1],[min(rog) max(rog)],'--k')
loglog(1000.*get_rvir(200).*[1 1],[min(rog) max(rog)],'--k')
loglog(2.*1000.*get_rvir(500).*[1 1],[min(rog) max(rog)],'--k')
xlabelmine('$r\,[\mathrm{kpc}]$');
ylabelmine('$\rho_{\mathrm{gas}}\,[\mathrm{M_{\odot}\,Mpc^{-3}}]$');
grid 
titlemine(sprintf('CL%s',num2str(cl)))

figure
loglog(rr.*1000,tp);
hold on 
loglog(1000.*get_rvir(500).*[1 1],[min(tp) max(tp)],'--k')
loglog(1000.*get_rvir(200).*[1 1],[min(tp) max(tp)],'--k')
loglog(2.*1000.*get_rvir(500).*[1 1],[min(tp) max(tp)],'--k')

xlabelmine('$r\,[\mathrm{kpc}]$');
ylabelmine('$T \,[\mathrm{K}]$');
grid
titlemine(sprintf('CL%s',num2str(cl)))

figure
loglog(rr.*1000,rog.*tp);
hold on 
loglog(1000.*get_rvir(500).*[1 1],[min(tp.*rog) max(tp.*rog)],'--k')
loglog(1000.*get_rvir(200).*[1 1],[min(tp.*rog) max(tp.*rog)],'--k')
loglog(2.*1000.*get_rvir(500).*[1 1],[min(tp.*rog) max(tp.*rog)],'--k')
titlemine(sprintf('CL%s',num2str(cl)))
xlabelmine('$r\,[\mathrm{kpc}]$');
ylabelmine('$P$');
grid 

pause


end