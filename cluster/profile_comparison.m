list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

plist=[1 0 1 0 1 1 1 0 0 1 0 0 0 1 1 0]; 
nplist=~plist;
plist=~nplist;

r=0.05:0.01:5; 

rop=zeros(length(list),length(r));
tp=zeros(length(list),length(r));
sp=zeros(length(list),length(r));

for i=1:length(list)
    new_env(sprintf('CL%d',list(i)),'csf','a1','win');
    [rop(i,:),~]=read_RHO_Profiles(r.*get_rvir);
    sp(i,:)=read_S_Profile(r.*get_rvir);
    tp(i,:)=read_T_Profile(r.*get_rvir);
end
ro_p=mean(rop(plist,:));
ro_np=mean(rop(nplist,:));
stro_p=std(rop(plist,:));
stro_np=std(rop(nplist,:));


s_p=mean(sp(plist,:));
s_np=mean(sp(nplist,:));
sts_p=std(sp(plist,:));
sts_np=std(sp(nplist,:));

global DEFAULT_PRINTOUT_DIR

t_p=mean(tp(plist,:));
t_np=mean(tp(nplist,:));
stt_p=std(tp(plist,:));
stt_np=std(tp(nplist,:));

figure
loglog(r,ro_p,'-b',r,ro_np,'-r','linewidth',1.5);
hold on
loglog(r,ro_p+0.5*stro_p,'--b',r,ro_np+0.5*stro_np,'--r');
loglog(r,ro_p-0.5*stro_p,'--b',r,ro_np-0.5*stro_np,'--r');
hold off
xlim([0.05 5])

hl=legend('Pen.','No Pen.');
set(hl,'Interpreter','Latex','Fontsize',12)
xlabelmine('$r/R_{vir}$');
ylabelmine('$\rho\,[\mathrm{M_\odot/Mpc^3}]$');
exportfig(gcf,sprintf('%s/density_pen_vs_nopen.png',DEFAULT_PRINTOUT_DIR),'FORMAT','png');


figure
loglog(r,t_p,'-b',r,t_np,'-r','linewidth',1.5);
legend('Pen.','No Pen.');
set(hl,'Interpreter','Latex','Fontsize',12)
xlabelmine('$r/R_{vir}$');
ylabelmine('$T\,[\mathrm{K}]$');
hold on
xlim([0.05 5])
loglog(r,t_p+0.5*stt_p,'--b',r,t_np+0.5*stt_np,'--r');
loglog(r,t_p-0.5*stt_p,'--b',r,t_np-0.5*stt_np,'--r');
hold off

exportfig(gcf,sprintf('%s/temperature_pen_vs_nopen.png',DEFAULT_PRINTOUT_DIR),'FORMAT','png');


figure
loglog(r,s_p,'-b',r,s_np,'-r','linewidth',1.5);
legend('Pen.','No Pen.');
set(hl,'Interpreter','Latex','Fontsize',12)
hold on
loglog(r,s_p+0.5*sts_p,'--b',r,s_np+0.5*sts_np,'--r');
loglog(r,s_p-0.5*sts_p,'--b',r,s_np-0.5*sts_np,'--r');
hold off
xlim([0.05 5])


xlabelmine('$r/R_{vir}$');
ylabelmine('$S\,[\mathrm{Kev\, cm^2}]$');


exportfig(gcf,sprintf('%s/entropy_pen_vs_nopen.png',DEFAULT_PRINTOUT_DIR),'FORMAT','png');
