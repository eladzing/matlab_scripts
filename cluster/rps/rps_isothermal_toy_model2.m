 output_dir='C:/Users/owner/Ubuntu One/cluster/printout';

%% basic plot, alpha=1   

 mu=[1e-5 1e-4 1e-3 1e-2 1e-1]; % virial mass ratio
rp=0.05:0.05:5; %in units of cluser rv
rs=zeros(length(mu),length(rp));
for i=1:length(mu)
rs(i,:)=mu(i).^(1/3).*rp./sqrt(1);
end

figure
plot(rp,rs,'Linewidth',2)

grid
ylim([1e-2 0.5])
hl=legend('$m_{sat}=10^{9}\,\mathrm{M_\odot}$','$m_{sat}=10^{10}\,\mathrm{M_\odot}$','$m_{sat}=10^{11}\,\mathrm{M_\odot}$','$m_{sat}=10^{12}\,\mathrm{M_\odot}$');%,'$m_{sat}=10^{13}\,\mathrm{M_\odot}$');
set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);
set(gca,'Fontsize',12);
xlabelmine('$r_p/R_c$')
ylabelmine('$r_s/R_s$')
titlemine('Stripping radius in a $10^{14}\,\mathrm{M_\odot}$ Cluster ($\alpha=1$)')
name='%s/rps_toy_m14.%s';
%exportfig(gcf,sprintf(name,output_dir,'png'),'format','png');
%exportfig(gcf,sprintf(name,output_dir,'eps'));


%% plot for exploring mu and alpha dependance 
figure 
mu=-5:0.01:-1.5; % virial mass ratio
alpha=0.01:0.01:1.0;

%rp=0.1:0.1:5; %in units of cluser rv
rp=1.0;
rs=zeros(length(mu),length(alpha));
for i=1:length(mu)
    for j=1:length(alpha)
        rs(i,j)=(10.^mu(i)).^(1/3).*rp./sqrt(alpha(j));
    end
end

map = brewermap(256,'RdBu');
colormap(map)
imagesc(mu,alpha,(1-(rs)').*100,[0 100]);
set(gca,'Ydir','normal','FontSize',14,'Ytick',0:0.1:1,'Xtick',-5:0.5:-1.5);
bar=colorbar;

barTitle(bar,'$M_{stripped}\,[\%]$');
%set(get(bar,'Title'),'String',,'Fontsize',12,'Interpreter','latex');
set(bar,'FontSize',12,'Ytick',0:10:100)
xlabelmine('$\log(\mu)$')
ylabelmine('$\varepsilon^{\prime}$')
%titlemine('Mass strippes at $r_p=1\,R_{\mathrm{vir}}$')

name='%s/rps_toy_mualf.%s';
%exportfig(gcf,sprintf(name,output_dir,'png'),'format','png');
%exportfig(gcf,sprintf(name,output_dir,'eps'));


% 
% %% plot fot alpha=1 
% figure
% plot(rp,rs,'Linewidth',2)
% 
% grid
% ylim([1e-2 0.5])
% hl=legend('$m_{sat}=10^{9}\,\mathrm{M_\odot}$','$m_{sat}=10^{10}\,\mathrm{M_\odot}$','$m_{sat}=10^{11}\,\mathrm{M_\odot}$','$m_{sat}=10^{12}\,\mathrm{M_\odot}$');%,'$m_{sat}=10^{13}\,\mathrm{M_\odot}$');
% set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);
% set(gca,'Fontsize',12);
% 
% titlemine('Stripping radius in a $10^{14}\,\mathrm{M_\odot}$ Cluster ($\alpha=1$)')
% %name='~/Ubuntu One/cluster/printout/rps_toy_model_14.%s';
% % exportfig(gcf,sprintf(name,'png'),'format','png');
% % exportfig(gcf,sprintf(name,'eps'));
% 
% %% plot fot alpha=0.5
% figure
% plot(rp,rs./sqrt(0.5),'Linewidth',2)
% 
% grid
% ylim([1e-2 0.5])
% hl=legend('$m_{sat}=10^{9}\,\mathrm{M_\odot}$','$m_{sat}=10^{10}\,\mathrm{M_\odot}$','$m_{sat}=10^{11}\,\mathrm{M_\odot}$','$m_{sat}=10^{12}\,\mathrm{M_\odot}$');%,'$m_{sat}=10^{13}\,\mathrm{M_\odot}$');
% set(hl,'Interpreter','latex','Location','NorthEast','Fontsize',14);
% set(gca,'Fontsize',12);
% xlabelmine('$r_p/R_c$')
% ylabelmine('$r_s/R_s$')
% titlemine('Stripping radius in a $10^{14}\,\mathrm{M_\odot}$ Cluster ($\alpha=0.5$)')
% %name='~/Ubuntu One/cluster/printout/rps_toy_model_14.%s';
% % exportfig(gcf,sprintf(name,'png'),'format','png');
% % exportfig(gcf,sprintf(name,'eps'));
% 
% %% plot fot alpha=0.5
% figure
% plot(rp,rs./sqrt(0.1),'Linewidth',2)
% 
% grid
% ylim([1e-2 0.5])
% hl=legend('$m_{sat}=10^{9}\,\mathrm{M_\odot}$','$m_{sat}=10^{10}\,\mathrm{M_\odot}$','$m_{sat}=10^{11}\,\mathrm{M_\odot}$','$m_{sat}=10^{12}\,\mathrm{M_\odot}$');%,'$m_{sat}=10^{13}\,\mathrm{M_\odot}$');
% set(hl,'Interpreter','latex','Location','Northeast','Fontsize',14);
% set(gca,'Fontsize',12);
% xlabelmine('$r_p/R_c$')
% ylabelmine('$r_s/R_s$')
% titlemine('Stripping radius in a $10^{14}\,\mathrm{M_\odot}$ Cluster ($\alpha=0.1$)')
% %name='~/Ubuntu One/cluster/printout/rps_toy_model_14.%s';
% % exportfig(gcf,sprintf(name,'png'),'format','png');
% % exportfig(gcf,sprintf(name,'eps'));
% 
% %% plot of et_S vs mass ration for rp=1,2,3 
% alpha=0.01:0.01:1;
% mu=[1e-5 1e-4 1e-3 1e-2]; % virial mass ratio
% %rp=0.01:0.01:5; %in units of cluser rv
% rs=zeros(length(mu),length(alpha));
% for i=1:length(mu)
% rs(i,:)=mu(i).^(1/3).*(1.0)./sqrt(2.*alpha);
% end
% 
% 
% 
% figure
% semilogx(alpha,rs','-','Linewidth',2)
% 
% grid
% ylim([0 1.0])
% %xlim([0.09 3])
% hl=legend('$m_{sat}=10^{9}\,\mathrm{M_\odot}$','$m_{sat}=10^{10}\,\mathrm{M_\odot}$','$m_{sat}=10^{11}\,\mathrm{M_\odot}$','$m_{sat}=10^{12}\,\mathrm{M_\odot}$');%,'$m_{sat}=10^{13}\,\mathrm{M_\odot}$');
% %hl=legend('$r_p=R_{vir}$','$r_p=2R_{vir}$');
% set(hl,'Interpreter','latex','Location','NorthEast','Fontsize',14);
% set(gca,'Fontsize',12);
% xlabelmine('$\alpha$')
% ylabelmine('$r_s/R_s$')
% titlemine('Stripping radius in a $10^{14}\,\mathrm{M_\odot}$ Cluster at $1\,R_{vir}$')
% %name='~/Ubuntu One/cluster/printout/rps_toy_model_14.%s';
% % exportfig(gcf,sprintf(name,'png'),'format','png');
% % exportfig(gcf,sprintf(name,'eps'));
