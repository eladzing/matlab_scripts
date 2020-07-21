%% this script uses the profiles constructed in batch_mk_profiles
% to create mean profiles of all halos.
% In addition the profiles are plotted

load('C:\Users\eladzing\OneDrive\cluster\matlab\mat_files\profiles.mat')

flist=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%mask=[ 0   0   1   0   0   1   0   0  0  1  0  0  0  1  1  0];
mask=true(size(flist));
%mask(1:2)=true;

rp=0.01:0.001:2.6; % common profile - in units of Rvir
new_env(101,'csf')
units;


% for jj=1:1
%    switch jj
%         case 1
%             lst=1:1:16;
%             printag='all';
%             titletag='All';
%                     case 2;
%                     lst=cflist(2,:);
%                     printag='cf'
%                     titletag='CF';
%                     case 3
%                     lst=ncflist(2,:);
%                     printag='ncf'
%                     titletag='NCF';
%     end

%% k=2   - for z=0.6  and   k=1   - for z=0
k=1;
nClst=length(profile(k).cluster);

n=0;
mdf=zeros(nClst,length(rp));
tp=zeros(nClst,length(rp));
rop=zeros(nClst,length(rp));
%vrp=zeros(nClst,length(rp));
sp=zeros(nClst,length(rp));

rhonorm=rho_mean(profile(k).zred).*deltavir(profile(k).zred);
for i=1:nClst
    %clname=sprintf('CL%d',cflist(i));
    prf=profile(k).cluster(i);
    n=n+1;
    
    mdf(i,:)=interp1(prf.rProflux,prf.fluxProf,rp); % units of Gyr^-1
    tp(i,:)=interp1(prf.rProf,prf.tmpProf,rp);
    rop(i,:)=interp1(prf.rProf,prf.rhoProf,rp)./rhonorm;
    %vrp(i,:)=interp1(stack{i,2},stack{i,6},rp,'linear');
    sp(i,:)=interp1(prf.rProf,prf.sProf,rp)./(kb/(1000*ev));
    
end

global DEFAULT_PRINTOUT_DIR
% printoutdir=

%%plot stacks
%figure;
xl=[2e-2 3];


subplot(4,2,1);
semilogx(rp,mdf','linewidth',1);grid;xlim(xl);ylim([-0.38 0.2]);
line([0.2 0.2],[-1e5 1e5],'Color',[0.2 0.8 0.2],'LineWidth',2);
%xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
ylabelmine('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$',12)
%title(sprintf('%s Clusters Mass Flux',titletag),'Fontsize',12);
%line([1e-3 1e1],[0 0],'Color','k');d

%legend(clname,'location','bestoutside');
set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'xticklabel','','Fontsize',12,'linewidth',1,...
    'Position',[0.15 0.77 0.42 0.21],'Fontsize',12);

fl1=mean(mdf,1);%fl2=sum(mdf2,1);
v1=std(mdf,0,1);%v2=std(mdf2,0,1);

subplot(4,2,2);
semilogx(rp,fl1,'-b','linewidth',2);grid;xlim(xl);ylim([-0.38 0.2]);hold on;
semilogx(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
hold off;

% xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
% ylabel('$\dot{M}/M_{\mathrm{vir}} [\mathrm{Gyr}^{-1}]$','Fontsize',12);

line([1e-3 1e1],[0 0],'Color','k');
line([0.2 0.2],[-1e5 1e5],'Color',[0.2 0.8 0.2],'LineWidth',2);
set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.4 -0.2 0 0.2],'xticklabel','','yticklabel','','Fontsize',12,'linewidth',1,...
    'Position',[0.57 0.77 0.42 0.21]);

%figure;
subplot(4,2,3);
loglog(rp,tp','linewidth',1);grid;xlim(xl);ylim([5e-2 10]);
line([0.2 0.2],[1e-5 1e5],'Color',[0.2 0.8 0.2],'linewidth',2);
%xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
ylabelmine('$T/T_{\mathrm{vir}}$');
set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'xticklabel','','Fontsize',12,'linewidth',1,...
    'Position',[0.15 0.56 0.42 0.21]);

fl1=mean(tp,1);%fl2=sum(mdf2,1);
v1=std(tp,0,1);%v2=std(mdf2,0,1);

%figure;
subplot(4,2,4);
loglog(rp,fl1,'-b','linewidth',2);grid;xlim(xl);ylim([5e-2 10]);hold on;
loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
line([0.2 0.2],[1e-5 1e5],'Color',[0.2 0.8 0.2],'linewidth',2);
hold off;
%xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
%ylabel('$T/T_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'xticklabel','','yticklabel','','Fontsize',12,'linewidth',1,...
    'Position',[0.57 0.56 0.42 0.21]);

%figure;
subplot(4,2,5)
loglog(rp,rop','linewidth',1);
grid;
xlim(xl);
ylim([1e-3 1e3]);
line([0.2 0.2],[1e-5 1e5],'Color',[0.2 0.8 0.2],'linewidth',2);
%xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
ylabelmine('$\rho_{gas}/\rho_{\mathrm{mean}}$')
set(gca,'XTick', [1e-2 1e-1 1 5],'YTick',[1e-2 1 1e2 ],'xticklabel','','Fontsize',12,'linewidth',1,...
    'Position',[0.15 0.35 0.42 0.21] );

fl1=mean(rop,1);%fl2=sum(mdf2,1);
v1=std(rop,0,1);%v2=std(mdf2,0,1);

%figure;
subplot(4,2,6);
loglog(rp,fl1,'-b','linewidth',2);grid;xlim(xl);ylim([1e-3 1e3]);hold on;
loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
line([0.2 0.2],[1e-5 1e5],'Color',[0.2 0.8 0.2],'linewidth',2);
hold off;
%xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
%ylabel('$\rho_{gas}/\rho_{\mathrm{vir}}$','Fontsize',12);
set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-2 1 1e2],'xticklabel','','yticklabel','','Fontsize',12,'linewidth',1,...
    'Position',[0.57 0.35 0.42 0.21]);

%figure;
subplot(4,2,7);
loglog(rp,sp','linewidth',1);grid;xlim(xl);ylim([1e1 1e4]);
line([0.2 0.2],[1e-5 1e5],'Color',[0.2 0.8 0.2],'linewidth',2);
xlabelmine('$r/R_{\mathrm{vir}}$')
ylabelmine('$S/k_B T_{vir} [\mathrm{cm^2}]$')
%set(gca,'XTick', [1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'Fontsize',12,'linewidth',1,...
set(gca,'XTick', [1e-2 1e-1 1 5],'Fontsize',12,'linewidth',1,...
'Position',[0.15 0.14 0.42 0.21]);


fl1=mean(sp,1);%fl2=sum(mdf2,1);
v1=std(sp,0,1);%v2=std(mdf2,0,1);

%figure;
subplot(4,2,8);
loglog(rp,fl1,'-b','linewidth',2);grid;xlim(xl);ylim([1e1 1e4]);hold on;
loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
line([0.2 0.2],[1e-5 1e5],'Color',[0.2 0.8 0.2],'linewidth',2);
hold off;
xlabelmine('$r/R_{\mathrm{vir}}$')
%ylabel('$S/S_{\mathrm{vir}}$','Fontsize',12);

set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'yticklabel','','Fontsize',12,'linewidth',1,...
    'Position',[0.57 0.14 0.42 0.21]);
%laprint(gcf,'ent_stk','options','factory','width',9,'factor',1,'scalefonts','on');

% if strcmp(pflag,'print')
%     %saveas(gcf,sprintf('%s/%s_stack_profs',result_dir,'multi'),'eps');
%     saveas(gcf,sprintf('%s/%s_stack_profs.png',result_dir,'multi'));
% end

%clear normv clname mdf tp rop vrp sp;
%end









