load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\profiles.mat')
rp=10.^(-2:0.001:log10(2.6));
new_env(101,'csf')
units;

k=1; % for z=0
nClst=length(profile(k).cluster);

n=0;
mdf=zeros(nClst,length(rp));


for i=1:nClst
    %clname=sprintf('CL%d',cflist(i));
    prf=profile(k).cluster(i);
    n=n+1;
    
    fac1=1;
    fac2=1; %0.11.*(prf.mvir/1e15).^0.15;
    %fac1/fac2
    
    mdf(i,:)=interp1(prf.rProflux,prf.fluxProf,rp).*fac1./fac2; % units of Gyr^-1
end

xl=[2e-2 3];

fl1=mean(mdf,1);%fl2=sum(mdf2,1);
fl2=median(mdf,1);
v1=std(mdf,0,1);%v2=std(mdf2,0,1);
fl1s=smooth(fl1,80);
uls=smooth(fl1+v1./2,80);
lls=smooth(fl1-v1./2,80);
%arY=cat(2,ul,fliplr(ll));
%arX=cat(2,rp,fliplr(rp));

figure
semilogx(rp,mdf(1:5,:)','linewidth',1.5);grid;xlim(xl);ylim([-0.3 0.15]);
hl=legend('CL101','CL102','CL103','CL104','CL105');
set(hl,'Interpreter','latex','fontsize',14,'Location','SouthWest')
%line([0.2 0.2],[-1e5 1e5],'Color',[0.2 0.8 0.2],'LineWidth',2);
set(gca,'fontsize',14)
xlabelmine('$r/R_{\mathrm{vir}}$',14)
ylabelmine('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$',14)

figure
semilogx(rp,mdf(6:10,:)','linewidth',1.5);grid;xlim(xl);ylim([-0.3 0.15]);
%line([0.2 0.2],[-1e5 1e5],'Color',[0.2 0.8 0.2],'LineWidth',2);
%xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
hl=legend('CL106','CL107','CL3','CL5','CL6');
set(hl,'Interpreter','latex','fontsize',14,'Location','SouthWest')

ylabelmine('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$',14)
xlabelmine('$r/R_{\mathrm{vir}}$',14)
set(gca,'fontsize',14)
printout_fig(gcf,'inflow2')


figure
semilogx(rp,mdf(11:16,:)','linewidth',1.5);grid;xlim(xl);ylim([-0.3 0.15]);
hl=legend('CL7','CL9','CL10','CL11','CL14','CL24');
set(hl,'Interpreter','latex','fontsize',14,'Location','SouthWest')

%line([0.2 0.2],[-1e5 1e5],'Color',[0.2 0.8 0.2],'LineWidth',2);
%xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
ylabelmine('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$',14)
xlabelmine('$r/R_{\mathrm{vir}}$',14)
set(gca,'fontsize',14)
printout_fig(gcf,'inflow3')

figure
semilogx(rp,fl1s,'-b','linewidth',2);grid;xlim(xl);ylim([-0.3 0.15]);hold on;
semilogx(rp,uls,'--b','linewidth',2);
semilogx(rp,lls,'--b','linewidth',2);
%semilogx(rp,fl2,':r','linewidth',2)
%semilogx(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-0.3 0.15]);
set(gca,'fontsize',14)

%line([0.2 0.2],[-1e5 1e5],'Color',[0.2 0.8 0.2],'LineWidth',2);
%xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
ylabelmine('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$',14)
xlabelmine('$r/R_{\mathrm{vir}}$',14)
printout_fig(gcf,'inflow_mean')

% 
% 
% %title(sprintf('%s Clusters Mass Flux',titletag),'Fontsize',12);
% %line([1e-3 1e1],[0 0],'Color','k');d
% 
% %legend(clname,'location','bestoutside');
% %set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'xticklabel','','Fontsize',12,'linewidth',1,...
% %    'Position',[0.15 0.77 0.42 0.21],'Fontsize',12);
% 
% fl1=mean(mdf,1);%fl2=sum(mdf2,1);
% v1=std(mdf,0,1);%v2=std(mdf2,0,1);

%subplot(4,2,2);
hold off;

% % xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
% % ylabel('$\dot{M}/M_{\mathrm{vir}} [\mathrm{Gyr}^{-1}]$','Fontsize',12);
% 
% line([1e-3 1e1],[0 0],'Color','k');
% line([0.2 0.2],[-1e5 1e5],'Color',[0.2 0.8 0.2],'LineWidth',2);
% %set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.4 -0.2 0 0.2],'xticklabel','','yticklabel','','Fontsize',12,'linewidth',1,...
% %    'Position',[0.57 0.77 0.42 0.21]);
