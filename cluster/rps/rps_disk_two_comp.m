%% plot for elucidating the force per unit area of gas in a disk for different
% disk types


beta=[0.2 0.5 1 2.0 5.0];  %ratio of Rg/Rd - extent of gas

r=0.001:0.001:20; % radial distance in units of Rd

gt=zeros(length(r),length(beta));
gs=zeros(length(r),length(beta));
gg=zeros(length(r),length(beta));

fg=0.05;
for i=1:length(beta)
    [gt(:,i) gs(:,i) gg(:,i)]=disk_accel_full(r,'beta',beta(i),'fg',fg);
    
    %frd(:,i)=disk_force_reduced(r,'beta',beta(i));
    %mf(:,i)=exp_disk_mass(r,beta(i));
end

figure
loglog(r,gt(:,1),r,gt(:,2),r,gt(:,3),r,gt(:,4),r,gt(:,5))

hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
set(hl,'Interpreter','latex','fontsize',12);
ylim([0.01 10])
xlim([0.01 20])
ylabelmine('$g(R)/\left(\pi G {\Sigma_0}\right)$')
xlabelmine('$R/R_d$')

figure
loglog(r,gg(:,1)./gt(:,1),r,gg(:,2)./gt(:,2),r,gg(:,3)./gt(:,3),r,gg(:,4)./gt(:,4),r,gg(:,5)./gt(:,5))

hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
set(hl,'Interpreter','latex','fontsize',12);
ylim([0.01 1])
xlim([0.01 20])
ylabelmine('$g_{gas}/g_{tot}$')
xlabelmine('$R/R_d$')






%
% if ~exist('plotflag','var')
%     plotflag=true(1,3);
% end

% %% plot mf vs fr
% if(plotflag(1))
%     figure
%     loglog(mf(:,1),frd(:,1),mf(:,2),frd(:,2),mf(:,3),frd(:,3),mf(:,4),frd(:,4),mf(:,5),frd(:,5));
%     
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     ylim([0.001 10])
%     xlim([0.01 1])
%     ylabelmine('$f(R)/\left(\pi G {\Sigma_0}^2 f_g\right)$')
%     xlabelmine('$M(R)/M_g$')
%     
%     figure
%     semilogx(mf(:,1),frd(:,1),mf(:,2),frd(:,2),mf(:,3),frd(:,3),mf(:,4),frd(:,4),mf(:,5),frd(:,5));
%     
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     xlim([0.01 1])
%     ylabelmine('$f(R)/\left(\pi G {\Sigma_0}^2 f_g\right)$')
%     xlabelmine('$M(R)/M_g$')
%     
%     
%     figure
%     plot(mf(:,1),frd(:,1),mf(:,2),frd(:,2),mf(:,3),frd(:,3),mf(:,4),frd(:,4),mf(:,5),frd(:,5));
%     xlim([0.01 1])
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     
%     ylabelmine('$f(R)/\left(\pi G {\Sigma_0}^2 f_g\right)$')
%     xlabelmine('$M(R)/M_g$')
% end
% 
% 
% 
% %% plot frd
% 
% if(plotflag(2))
%     figure
%     loglog(r,frd(:,1),r,frd(:,2),r,frd(:,3),r,frd(:,4),r,frd(:,5));
%     
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     ylim([0.001 1.1])
%     xlabelmine('$R/R_d$')
%     ylabelmine('$f(R)/\left(\pi G {\Sigma_0}^2 f_g\right)$')
%     
%     figure
%     semilogx(r,frd(:,1),r,frd(:,2),r,frd(:,3),r,frd(:,4),r,frd(:,5));
%     
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     xlim([0 5])
%     xlabelmine('$R/R_d$')
%     ylabelmine('$f(R)/\left(\pi G {\Sigma_0}^2 f_g\right)$')
%     
%     
%     figure
%     plot(r,frd(:,1),r,frd(:,2),r,frd(:,3),r,frd(:,4),r,frd(:,5));
%     
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     xlim([0 5])
%     xlabelmine('$R/R_d$')
%     ylabelmine('$f(R)/\left(\pi G {\Sigma_0}^2 f_g\right)$')
% end
% 
% 
% %% plot mass profiles
% if(plotflag(3))
%     figure
%     loglog(r,mf(:,1),r,mf(:,2),r,mf(:,3),r,mf(:,4),r,mf(:,5));
%     
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     ylim([0.001 10])
%     xlabelmine('$R/R_d$')
%     ylabelmine('$M(R)/M_g$')
%     
%     
%     figure
%     semilogx(r,mf(:,1),r,mf(:,2),r,mf(:,3),r,mf(:,4),r,mf(:,5));
%     
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     xlim([0.01 1])
%     xlabelmine('$R/R_d$')
%     ylabelmine('$M(R)/M_g$')
%     
%     
%     figure
%     plot(r,mf(:,1),r,mf(:,2),r,mf(:,3),r,mf(:,4),r,mf(:,5));
%     
%     hl=legend('$\beta=0.2$','$\beta=0.5$','$\beta=1$','$\beta=2$','$\beta=5$');
%     set(hl,'Interpreter','latex','fontsize',12);
%     xlim([0.01 1])
%     xlabelmine('$R/R_d$')
%     ylabelmine('$M(R)/M_g$')
% end