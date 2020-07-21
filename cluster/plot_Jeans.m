%% plot stacked jeans

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
load('C:\Users\eladzing\OneDrive\cluster\matlab\mat_files\JeansShell.mat')

aexp='a06';
penetration_depth_a1


switch aexp
    case 'a1'
        jean=Jean_a1;
        rlxMask=rlx;
        urlxMask=urlx;
    case 'a06'
        jean=Jean_a06;
        rlxMask=rlx06;
        urlxMask=urlx06;
    otherwise
        error('PLOT_JEANS: wrong aexp %s',aexp)
end
        


%mask=~(maxpen06<=25);
%mask=true(size(list));

r0=logspace(-2,2,2000);%       0.01:0.001:1;%     logspace(-2,0.3,1000);%

therm=zeros(length(list),length(r0));
turb=zeros(length(list),length(r0));
spher=zeros(length(list),length(r0));
isotr=zeros(length(list),length(r0));
vPress=zeros(length(list),length(r0));
tot=zeros(length(list),length(r0));

for i=1:length(jean)
    
    
    
    gg0=jean(i).grav;
    th0=jean(i).fTherm;
    tu0=jean(i).fTurb;
    iso0=jean(i).isotrop(2:end-1);
    sph0=jean(i).spheric(2:end-1);
    
    rp=jean(i).rp(2:end-1)./jean(i).rv;
    
    gg=interp1(rp,gg0,r0);
    th=interp1(rp,th0,r0);
    tu=interp1(rp,tu0,r0);
    iso=interp1(rp,iso0,r0);
    sph=interp1(rp,sph0,r0);
    
    
    
    therm(i,:)=(th./gg);
    turb(i,:)=(tu./gg);
    spher(i,:)=-(sph./gg);
    isotr(i,:)=(iso./gg);
    vPress(i,:)=((tu+iso-sph)./gg);
    tot(i,:)=(th+tu+iso-sph)./gg;
    %tot(i,:)=(th+tu)./gg;
    
     %% plot individual
     hf=figure;
     h=[];
     h(1)=semilogx(r0,th./gg,'-r','DisplayName','$T/g$')  ;
     hold on
     h(2)=semilogx(r0,tu./gg,'-b','DisplayName','$\sigma/g$')  ;
     h(3)=semilogx(r0,iso./gg,'-g','DisplayName','$iso/g$')  ;
     h(4)=semilogx(r0,-sph./gg,'-c','DisplayName','$sph$')  ;
     h(5)=semilogx(r0,(th+tu+iso-sph)./gg,'-k','DisplayName','$Tot$')  ;
     semilogx([0.01 2],[1 1],':k')
     hl=clickableLegend(h);
     set(hl,'Interpreter','latex','fontsize',14,'Location','NorthWest')
     xlim([0.01 2])
     ylim([0 3])
     xlabelmine('$r/R_{\mathrm{vir}}$')
     titlemine(sprintf('CL%s',num2str(jean(i).cl)));
     %pause
    
end

%% sum up profile
xl=[0.03 2.0];
yl1=[0 1.6];
yl2=[-0.2 0.8];
nSmooth=200;
meth='moving';

%% rlx 
mask=rlxMask;
thrm=smooth(r0,mean(therm(mask,:),1),nSmooth,meth);
nonThrm=smooth(r0,mean(vPress(mask,:),1),nSmooth,meth);
trb=smooth(r0,mean(turb(mask,:),1),nSmooth,meth);
iso=smooth(r0,mean(isotr(mask,:),1),nSmooth,meth);
sph=smooth(r0,mean(spher(mask,:),1),nSmooth,meth);
tol=smooth(r0,mean(tot(mask,:),1),nSmooth,meth);

figure
%subplot(1,12,1:5)
h(1)=semilogx(r0,thrm,'r','linewidth',2,'DisplayName','$F_{Th}/F_g$');
hold on
h(2)=semilogx(r0,nonThrm,'b','linewidth',2,'DisplayName','$F_{NT}/F_g$');
%h(3)=semilogx(r0,trb,'g','linewidth',2,'DisplayName','$P_{\sigma}/g$');
%h(4)=semilogx(r0,iso,'c','linewidth',2,'DisplayName','Isotrop.');
%h(5)=semilogx(r0,sph,'m','linewidth',2,'DisplayName','Spher.');
h(6)=semilogx(r0,tol,'k','linewidth',2,'DisplayName','$(F_{TH}+F_{NT})/F_g$');
semilogx(r0,ones(size(r0)),'--k','linewidth',2)
xlabelmine('$r/R_{\mathrm{vir}}$');
xlim(xl)
ylim(yl1)
grid
hl=legend(h([1 2 6]));
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthWest')
%titlemine(sprintf('Relaxed %s',aexp))
set(gca,'fontsize',14,'YTick',0:0.2:2)
printout_fig(gcf,sprintf('jeans_rlx_%s',aexp));

figure
%subplot(1,12,1:5)
h(2)=semilogx(r0,nonThrm,'b','linewidth',2,'DisplayName','$F_{NT}/F_g$');
hold on
h(3)=semilogx(r0,trb,'g','linewidth',2,'DisplayName','$F_{\sigma}/F_g$');
h(4)=semilogx(r0,iso,'c','linewidth',2,'DisplayName','Anisotropy');
h(5)=semilogx(r0,sph,'m','linewidth',2,'DisplayName','Asphericity');
%h(6)=semilogx(r0,tol,'k','linewidth',2,'DisplayName','Total');
semilogx(r0,ones(size(r0)),'--k','linewidth',2)
xlabelmine('$r/R_{\mathrm{vir}}$');
xlim(xl)
ylim(yl2)
grid
hl=legend(h(2:5));
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthWest')
%titlemine(sprintf('Relaxed %s',aexp))
set(gca,'fontsize',14,'YTick',-0.4:0.1:2)
printout_fig(gcf,sprintf('jeansNT_rlx_%s',aexp));



%% urlx 
mask=urlxMask;
thrm=smooth(r0,mean(therm(mask,:),1),nSmooth,meth);
nonThrm=smooth(r0,mean(vPress(mask,:),1),nSmooth,meth);
trb=smooth(r0,mean(turb(mask,:),1),nSmooth,meth);
iso=smooth(r0,mean(isotr(mask,:),1),nSmooth,meth);
sph=smooth(r0,mean(spher(mask,:),1),nSmooth,meth);
tol=smooth(r0,mean(tot(mask,:),1),nSmooth,meth);

figure
%subplot(1,12,1:5)
h(1)=semilogx(r0,thrm,'r','linewidth',2,'DisplayName','$F_{Th}/F_g$');
hold on
h(2)=semilogx(r0,nonThrm,'b','linewidth',2,'DisplayName','$F_{NT}/F_g$');
%h(3)=semilogx(r0,trb,'g','linewidth',2,'DisplayName','$P_{\sigma}/g$');
%h(4)=semilogx(r0,iso,'c','linewidth',2,'DisplayName','Isotrop.');
%h(5)=semilogx(r0,sph,'m','linewidth',2,'DisplayName','Spher.');
h(6)=semilogx(r0,tol,'k','linewidth',2,'DisplayName','$(F_{TH}+F_{NT})/F_g$');
semilogx(r0,ones(size(r0)),'--k','linewidth',2)
xlabelmine('$r/R_{\mathrm{vir}}$');
xlim(xl)
ylim(yl1)
grid
hl=legend(h([1 2 6]));
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthWest')
%titlemine(sprintf('Relaxed %s',aexp))
set(gca,'fontsize',14,'YTick',0:0.2:2)
printout_fig(gcf,sprintf('jeans_urlx_%s',aexp));

figure
%subplot(1,12,1:5)
h(2)=semilogx(r0,nonThrm,'b','linewidth',2,'DisplayName','$F_{NT}/F_g$');
hold on
h(3)=semilogx(r0,trb,'g','linewidth',2,'DisplayName','$F_{\sigma}/F_g$');
h(4)=semilogx(r0,iso,'c','linewidth',2,'DisplayName','Anisotropy');
h(5)=semilogx(r0,sph,'m','linewidth',2,'DisplayName','Asphericity');
%h(6)=semilogx(r0,tol,'k','linewidth',2,'DisplayName','Total');
semilogx(r0,ones(size(r0)),'--k','linewidth',2)
xlabelmine('$r/R_{\mathrm{vir}}$');
xlim(xl)
ylim(yl2)
grid
hl=legend(h(2:5));
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthWest')
%titlemine(sprintf('Relaxed %s',aexp))
set(gca,'fontsize',14,'YTick',-0.4:0.1:2)
printout_fig(gcf,sprintf('jeansNT_urlx_%s',aexp));



% 
% 
% 
% figure
% %subplot(1,12,6:12)
% 
% h(1)=semilogx(r0,thrm,'r','linewidth',2,'DisplayName','$P_{T}/g$');
% hold on
% h(2)=semilogx(r0,nonThrm,'b','linewidth',2,'DisplayName','$P_{NT}/g$');
% h(3)=semilogx(r0,trb,'g','linewidth',2,'DisplayName','$P_{\sigma}/g$');
% h(4)=semilogx(r0,iso,'c','linewidth',2,'DisplayName','Iso.');
% h(5)=semilogx(r0,sph,'m','linewidth',2,'DisplayName','Sph.');
% h(6)=semilogx(r0,tol,'k','linewidth',2,'DisplayName','Tot.');
% semilogx(r0,ones(size(r0)),'--k','linewidth',2)
% xlabelmine('$r/R_{\mathrm{vir}}$');
% xlim(xl)
% ylim([-0.2 1.6])
% grid
% hl=legend(h); %hl=gridLegend(h); %clickableLegend(h);
% set(hl,'Interpreter','latex','fontsize',12,'Location','BestOutside')
% %set(hl,'fontsize',14,'Location','Best')
% titlemine(sprintf('UnRelaxed %s',aexp))
% set(gca,'fontsize',14,'YTick',0:0.2:1.4)
% printout_fig(gcf,sprintf('jeans_urlx_%s',aexp));
