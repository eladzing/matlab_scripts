%% this script uses the profiles constructed in batch_mk_profiles
% to create mean profiles of all halos.
% In addition the profiles are plotted


%function plot_grouped_profiles(aexp)
load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\profiles.mat')
penetration_depth_a1
flist=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%mask=[ 0   0   1   0   0   1   0   0  0  1  0  0  0  1  1  0];
mask=true(size(flist));
%mask(1:2)=true;

%rp=0.01:0.001:2.6; % common profile - in units of Rvir
rp=10.^(-2:0.01:log10(2.6));
k=1;
switch(aexp)
    case 'a1'
        new_env(101,'csf','a1')
        k=1;
        rMask=rlx;
        uMask=urlx;
        pMask=maxpen<25;
        nMask=~pMask;
        ylF=[-0.2 0.2];
        ylEZ=[30 150];
    case 'a06'
        new_env(101,'csf','a06')
        k=2;
        rMask=rlx06;
        uMask=urlx06;
        pMask=maxpen06<25;
        nMask=~pMask;
        ylF=[-0.7 0.7];
        ylEZ=[15 70];
    otherwise
        error('plot_grouped_profiles: wrong aexp %s',aexp);
end
units;
aMask=true(size(flist));
%% k=2   - for z=0.6  and   k=1   - for z=0



global zred
nClst=length(flist);

n=0;
mdf=zeros(nClst,length(rp));
tp=zeros(nClst,length(rp));
rop=zeros(nClst,length(rp));
%vrp=zeros(nClst,length(rp));
sp=zeros(nClst,length(rp));

rhonorm=rho_mean(zred).*deltavir(zred);


for i=1:nClst
    %clname=sprintf('CL%d',cflist(i));
    prf=profileInOut(k).cluster(i);
    n=n+1;
    
    mdf(i,:)=interp1(prf.rProflux,prf.fluxProf,rp); % units of Gyr^-1
    mdfIn(i,:)=interp1(prf.rProflux,prf.fluxProfIn,rp); % units of Gyr^-1
    mdfOut(i,:)=interp1(prf.rProflux,prf.fluxProfOut,rp); % units of Gyr^-1
    
    tp(i,:)=interp1(prf.rProf,prf.tmpProf,rp);
    rop(i,:)=interp1(prf.rProf,prf.rhoProf,rp)./rhonorm;
    %vrp(i,:)=interp1(stack{i,2},stack{i,6},rp,'linear');
    sp(i,:)=interp1(prf.rProf,prf.sProf,rp)./(kb/(1000*ev));
    
end

%% plot 

%global DEFAULT_PRINTOUT_DIR
% printoutdir=

%%plot stacks
%figure;
xl=[3e-2 2.9];
xlZoom=[3e-2 2e-1];

%% flux stuff 
nSmooth=10;

lFluxIA=aux(mdfIn,aMask,nSmooth);
lFluxIR=aux(mdfIn,rMask,nSmooth);
lFluxIU=aux(mdfIn,uMask,nSmooth);
lFluxIP=aux(mdfIn,pMask,nSmooth);
lFluxIN=aux(mdfIn,nMask,nSmooth);

lFluxOA=aux(mdfOut,aMask,nSmooth);
lFluxOR=aux(mdfOut,rMask,nSmooth);
lFluxOU=aux(mdfOut,uMask,nSmooth);
lFluxOP=aux(mdfOut,pMask,nSmooth);
lFluxON=aux(mdfOut,nMask,nSmooth);


figure 
h=[];
h(1)=semilogx(rp,lFluxIR.meen,'-b','linewidth',2,'DisplayName','Relaxed In');
hold on
semilogx(rp,lFluxIR.top,'--b','linewidth',2,'DisplayName','Relaxed In');
semilogx(rp,lFluxIR.bot,'--b','linewidth',2,'DisplayName','Relaxed In');

h(2)=semilogx(rp,lFluxOR.meen,'-r','linewidth',2,'DisplayName','Relaxed Out');
semilogx(rp,lFluxOR.top,'--r','linewidth',2,'DisplayName','Relaxed In');
semilogx(rp,lFluxOR.bot,'--r','linewidth',2,'DisplayName','Relaxed In');
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthWest')

grid;xlim(xl);ylim(ylF);
set(gca,'fontsize',14)
xlabelmine('$r/R_{\mathrm{vir}}$')
ylabelmine('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$',14)

printout_fig(gcf,sprintf('stkProfs_flux_rlx_%s',aexp),'v');


figure
h=[];
h(1)=semilogx(rp,lFluxIU.meen,'-b','linewidth',2,'DisplayName','Unrelax In');
hold on 
semilogx(rp,lFluxIU.top,'--b','linewidth',2,'DisplayName','Unrelax');
semilogx(rp,lFluxIU.bot,'--b','linewidth',2,'DisplayName','Unrelax');

h(2)=semilogx(rp,lFluxOU.meen,'-r','linewidth',2,'DisplayName','Unrelax Out');
semilogx(rp,lFluxOU.top,'--r','linewidth',2,'DisplayName','Unrelax Out');
semilogx(rp,lFluxOU.bot,'--r','linewidth',2,'DisplayName','Unrelax Out');

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthWest')

grid;xlim(xl);ylim(ylF);
set(gca,'fontsize',14)
xlabelmine('$r/R_{\mathrm{vir}}$')
ylabelmine('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$',14)
printout_fig(gcf,sprintf('stkProfs_flux_urlx_%s',aexp),'v');


%% temperature 

nSmooth=10;

lTempA=aux(tp,aMask,nSmooth);
lTempR=aux(tp,rMask,nSmooth);
lTempU=aux(tp,uMask,nSmooth);
lTempP=aux(tp,pMask,nSmooth);
lTempN=aux(tp,nMask,nSmooth);


figure 
h=[];
h(1)=loglog(rp,lTempR.meen,'-b','linewidth',2,'DisplayName','Relaxed');
hold on
loglog(rp,lTempR.top,'--b','linewidth',2,'DisplayName','Relaxeded');
loglog(rp,lTempR.bot,'--b','linewidth',2,'DisplayName','Relaxeded');

h(2)=loglog(rp,lTempU.meen,'-r','linewidth',2,'DisplayName','Unrelaxed');
loglog(rp,lTempU.top,'--r','linewidth',2,'DisplayName','Unrelaxed');
loglog(rp,lTempU.bot,'--r','linewidth',2,'DisplayName','Unrelaxed');
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthEast')

grid;xlim(xl);ylim([0.1 3.1]);
grid minor
set(gca,'fontsize',14)
xlabelmine('$r/R_{\mathrm{vir}}$')
ylabelmine('$T/T_{\mathrm{vir}}$');

% insert zoom box 
axes2 = axes('ZMinorGrid','on','YScale','log',...
    'YMinorTick','on',...
    'YMinorGrid','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'XMinorGrid','on',...
    'Position',[0.2 0.2 0.3 0.37],...
    'FontSize',12);
box(axes2,'on');
grid(axes2,'on');
hold(axes2,'all');

loglog(rp,lTempR.meen,'-b','linewidth',2,'DisplayName','Relaxed');
loglog(rp,lTempR.top,'--b','linewidth',2,'DisplayName','Relaxeded');
loglog(rp,lTempR.bot,'--b','linewidth',2,'DisplayName','Relaxeded');
loglog(rp,lTempU.meen,'-r','linewidth',2,'DisplayName','Unrelaxed');
loglog(rp,lTempU.top,'--r','linewidth',2,'DisplayName','Unrelaxed');
loglog(rp,lTempU.bot,'--r','linewidth',2,'DisplayName','Unrelaxed');
xlim(axes2,[0.03 0.3]);
ylim(axes2,[1.0 3.1])

printout_fig(gcf,sprintf('stkProfs_temperature_%s',aexp),'v');
% % zoom in
% figure 
% h=[];
% h(1)=loglog(rp,lTempR.meen,'-b','linewidth',2,'DisplayName','Relaxed');
% hold on
% loglog(rp,lTempR.top,'--b','linewidth',2,'DisplayName','Relaxeded');
% loglog(rp,lTempR.bot,'--b','linewidth',2,'DisplayName','Relaxeded');
% 
% h(2)=loglog(rp,lTempU.meen,'-r','linewidth',2,'DisplayName','Unrelaxed');
% loglog(rp,lTempU.top,'--r','linewidth',2,'DisplayName','Unrelaxed');
% loglog(rp,lTempU.bot,'--r','linewidth',2,'DisplayName','Unrelaxed');
% hl=legend(h);
% set(hl,'Interpreter','latex','fontsize',14,'Location','NorthEast')
% 
% grid;xlim(xlZoom);ylim([1 5]);
% grid minor
% set(gca,'fontsize',14)
% xlabelmine('$r/R_{\mathrm{vir}}$')
% ylabelmine('$T/T_{\mathrm{vir}}$');


%% density 

nSmooth=10;

lRhoA=aux(rop,aMask,nSmooth);
lRhoR=aux(rop,rMask,nSmooth);
lRhoU=aux(rop,uMask,nSmooth);
lRhoP=aux(rop,pMask,nSmooth);
lRhoN=aux(rop,nMask,nSmooth);


figure 
h=[];
h(1)=loglog(rp,lRhoR.meen,'-b','linewidth',2,'DisplayName','Relaxed');
hold on
loglog(rp,lRhoR.top,'--b','linewidth',2,'DisplayName','Relaxeded');
loglog(rp,lRhoR.bot,'--b','linewidth',2,'DisplayName','Relaxeded');

h(2)=loglog(rp,lRhoU.meen,'-r','linewidth',2,'DisplayName','Unrelaxed');
loglog(rp,lRhoU.top,'--r','linewidth',2,'DisplayName','Unrelaxed');
loglog(rp,lRhoU.bot,'--r','linewidth',2,'DisplayName','Unrelaxed');
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthEast')

grid;xlim(xl);ylim([0.001 30]);
grid minor
set(gca,'fontsize',14)
xlabelmine('$r/R_{\mathrm{vir}}$')
ylabelmine('$\rho_{gas}/\rho_{\mathrm{vir}}$')


axes2 = axes('ZMinorGrid','on','YScale','log',...
    'YMinorTick','on',...
    'YMinorGrid','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'XMinorGrid','on',...
    'Position',[0.2 0.2 0.3 0.37],...
    'FontSize',12);
box(axes2,'on');
grid(axes2,'on');
hold(axes2,'all');

loglog(rp,lRhoR.meen,'-b','linewidth',2,'DisplayName','Relaxed');
loglog(rp,lRhoR.top,'--b','linewidth',2,'DisplayName','Relaxeded');
loglog(rp,lRhoR.bot,'--b','linewidth',2,'DisplayName','Relaxeded');
loglog(rp,lRhoU.meen,'-r','linewidth',2,'DisplayName','Unrelaxed');
loglog(rp,lRhoU.top,'--r','linewidth',2,'DisplayName','Unrelaxed');
loglog(rp,lRhoU.bot,'--r','linewidth',2,'DisplayName','Unrelaxed');
xlim(axes2,[0.03 0.3]);
ylim(axes2,[1.0 20])

printout_fig(gcf,sprintf('stkProfs_density_%s',aexp),'v');


