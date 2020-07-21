%% test dependence on delta vir

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

for i=1:length(list)
    new_env(list(i),'a06')
    mv(i)=get_mvir./1e14;
    m200(i)=get_mvir(200)./1e14;
    m500(i)=get_mvir(500)./1e14;
    rv(i)=get_rvir;
    r200(i)=get_rvir(200);
    r500(i)=get_rvir(500);
end
global zred
[mv ix]=sort(mv);
m200=m200(ix);
m500=m500(ix);
rv=rv(ix);
r200=r200(ix);
r500=r500(ix);




[fitR200 gofR200]=linearFit(rv,r200,'names','$R_{\mathrm{vir}}$','$R_{200}$','noshow');
[fitR500 gofR500]=linearFit(rv,r500,'names','$R_{\mathrm{vir}}$','$R_{500}$','noshow');
[fitM200 gofM200]=linearFit(mv,m200,'names','$M_{\mathrm{vir}}$','$M_{200}$','noshow');
[fitM500 gofM500]=linearFit(mv,m500,'names','$M_{\mathrm{vir}}$','$M_{500}$','noshow');

ratR200=r200./rv;
ratR500=r500./rv;
ratM200=m200./mv;
ratM500=m500./mv;

[fitRatR200 gofRatR200]=linearFit(mv,ratR200,'names','$M_{\mathrm{vir}}$','$R_{200}/R_{\mathrm{vir}}$','noshow');
[fitRatR500 gofRatR500]=linearFit(mv,ratR500,'names','$M_{\mathrm{vir}}$','$R_{500}/R_{\mathrm{vir}}$','noshow');
[fitRatM200 gofRatM200]=linearFit(mv,ratM200,'names','$M_{\mathrm{vir}}$','$M_{200}/M_{\mathrm{vir}}$','noshow');
[fitRatM500 gofRatM500]=linearFit(mv,ratM500,'names','$M_{\mathrm{vir}}$','$M_{500}/M_{\mathrm{vir}}$','noshow');

ll=1:length(list);
[fit2RatR200 gof2RatR200]=linearFit(ll,ratR200,'names','$R_{\mathrm{vir}}$','$R_{200}/R_{\mathrm{vir}}$','noshow');
[fit2RatR500 gof2RatR500]=linearFit(ll,ratR500,'names','$R_{\mathrm{vir}}$','$R_{500}/R_{\mathrm{vir}}$','noshow');
[fit2RatM200 gof2RatM200]=linearFit(ll,ratM200,'names','$R_{\mathrm{vir}}$','$M_{200}/M_{\mathrm{vir}}$','noshow');
[fit2RatM500 gof2RatM500]=linearFit(ll,ratM500,'names','$R_{\mathrm{vir}}$','$M_{500}/M_{\mathrm{vir}}$','noshow');



afR200=fitR200.p1;
afR500=fitR500.p1;
afM200=fitM200.p1;
afM500=fitM500.p1;

lRatR200=[fitRatR200.p1 fitRatR200.p2];
lRatR500=[fitRatR500.p1 fitRatR500.p2];
lRatM200=[fitRatM200.p1 fitRatM200.p2];
lRatM500=[fitRatM500.p1 fitRatM500.p2];

l2RatR200=[fit2RatR200.p1 fit2RatR200.p2];
l2RatR500=[fit2RatR500.p1 fit2RatR500.p2];
l2RatM200=[fit2RatM200.p1 fit2RatM200.p2];
l2RatM500=[fit2RatM500.p1 fit2RatM500.p2];


a200=1/sqrt(200/deltavir(zred));
a500=1/sqrt(500/deltavir(zred));

l200=rv.*a200;%+r200(rv==min(rv));
l500=rv.*a500; %+r500(rv==min(rv));
lf200=rv.*afR200;%+r200(rv==min(rv));
lf500=rv.*afR500; %+r500(rv==min(rv));



hf1=figure;
h=[];
h(1)=plot(rv,r200,'+b','DisplayName','$R_{200}$','MarkerSize',12,'Linewidth',2);
hold on
h(2)=plot(rv,r500,'+r','DisplayName','$R_{500}$','MarkerSize',12,'Linewidth',2);
h(3)=plot(rv,l200,'--b','Linewidth',2,'DisplayName','$\propto \sqrt{200/\Delta_{\mathrm{vir}}}$');
h(4)=plot(rv,l500,'--r','Linewidth',2,'DisplayName','$\propto \sqrt{500/\Delta_{\mathrm{vir}}}$');
%h(5)=plot(rv,lf200,'-b','Linewidth',2,'DisplayName','fit');
%h(6)=plot(rv,lf500,'-r','Linewidth',2,'DisplayName','fit');



plot([0 4],[0 4],'--k','linewidth',2)
hold off
xlabelmine('$R_{\mathrm{vir}}\,[\mathrm{Mpc}]$')
ylabelmine('$R_{200},R_{500}\,[\mathrm{Mpc}]$')
grid
hl=legend(h);
%set(hl,'Interpreter','latex','fontsize',12,'Location','SouthEast');
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthWest');

set(gca,'fontsize',14,'box','on')

%xlim([1 4])
%ylim([1 4])


xlim([0.4 1.8])
ylim([0.4 1.8])
axis square
printout_fig(gcf,'rvir_vs_r200_vs_r500_a06')



l200=mv.*a200;%+r200(rv==min(rv));
l500=mv.*a500; %+r500(rv==min(rv));
lf200=mv.*afM200;%+r200(rv==min(rv));
lf500=mv.*afM500; %+r500(rv==min(rv));


h=[];
h(1)=loglog(mv,m200,'+b','DisplayName','$M_{200}$','MarkerSize',10,'Linewidth',2);
hold on
h(2)=loglog(mv,m500,'+r','DisplayName','$M_{500}$','MarkerSize',10,'Linewidth',2);
h(3)=loglog(mv,l200,'--b','Linewidth',2,'DisplayName','$\propto \sqrt{200/\Delta_{\mathrm{vir}}}$');
h(4)=loglog(mv,l500,'--r','Linewidth',2,'DisplayName','$\propto \sqrt{500/\Delta_{\mathrm{vir}}}$');
h(5)=loglog(mv,lf200,'-b','Linewidth',2,'DisplayName','fit');
h(6)=loglog(mv,lf500,'-r','Linewidth',2,'DisplayName','fit');

xlim([0.6 30])
ylim([0.6 30])

loglog(xlim,ylim,':k')
hold off
xlabelmine('$M_{\mathrm{vir}}\,[10^{14}\mathrm{M_\odot}]$')
ylabelmine('$M_{200},M_{500}\,[\mathrm{M_\odot}]$')
grid
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthWest');

set(gca,'fontsize',12,'box','on')
%axis square

%printout_fig(gcf,'mvir_vs_m200_vs_m500')

hf3=figure;

h=[];
h(1)=semilogx(mv,ratR200,'ob','DisplayName','$R_{200}/R_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);
hold on
h(2)=semilogx(mv,ratR500,'or','DisplayName','$R_{500}/R_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);
h(3)=semilogx(mv,ratM200,'db','DisplayName','$M_{200}/M_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);
h(4)=semilogx(mv,ratM500,'dr','DisplayName','$M_{500}/M_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);

h(5)=semilogx(mv,a200.*ones(size(mv)),'--b','Linewidth',2,'DisplayName','$\sqrt{200/\Delta_{\mathrm{vir}}}$');
h(6)=semilogx(mv,a500.*ones(size(mv)),'--r','Linewidth',2,'DisplayName','$\sqrt{500/\Delta_{\mathrm{vir}}}$');

xlim([0.6 30])
%ylim([0.6 30])

% 

h(7)=semilogx(xlim,xlim.*lRatR200(1)+lRatR200(2),':b','Linewidth',2,'DisplayName','fit');
h(8)=semilogx(xlim,xlim.*lRatR500(1)+lRatR500(2),':r','Linewidth',2,'DisplayName','fit');
h(9)=semilogx(xlim,xlim.*lRatM200(1)+lRatM200(2),':b','Linewidth',2,'DisplayName','fit');
h(10)=semilogx(xlim,xlim.*lRatM500(1)+lRatM500(2),':r','Linewidth',2,'DisplayName','fit');


hold off
xlabelmine('$M_{\mathrm{vir}}\,\mathrm{10^{14}M_\odot}]$')
ylabelmine('Ratio');%            '$M_{200},M_{500}\,[\mathrm{M_\odot}]$')
grid
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEastOutside');

set(gca,'fontsize',12,'box','on')




hf4=figure;

h=[];
h(1)=semilogx(mv,ratR200,'ob','DisplayName','$R_{200}/R_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);
hold on
h(2)=semilogx(mv,ratM200,'dr','DisplayName','$M_{200}/M_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);


xlim([0.6 30])
ylim([1.02 1.32])

% 
h(3)=semilogx(xlim,a200.*ones(size(xlim)),'--k','Linewidth',2,'DisplayName','$\sqrt{200/\Delta_{\mathrm{vir}}}$');
h(4)=semilogx(xlim,xlim.*lRatR200(1)+lRatR200(2),'-b','Linewidth',2,'DisplayName','fit');
h(5)=semilogx(xlim,xlim.*lRatM200(1)+lRatM200(2),'-r','Linewidth',2,'DisplayName','fit');

hold off
xlabelmine('$M_{\mathrm{vir}}\,\mathrm{10^{14}M_\odot}]$')
ylabelmine('Ratio');%            '$M_{200},M_{500}\,[\mathrm{M_\odot}]$')
grid
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEastOutside');

set(gca,'fontsize',12,'box','on')


hf5=figure;

h=[];

h(1)=semilogx(mv,ratR500,'ob','DisplayName','$R_{500}/R_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);
hold on
h(2)=semilogx(mv,ratM500,'dr','DisplayName','$M_{500}/M_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);


xlim([0.6 30])
ylim([0.77 0.97])

% 
h(3)=semilogx(xlim,a500.*ones(size(xlim)),'--k','Linewidth',2,'DisplayName','$\sqrt{500/\Delta_{\mathrm{vir}}}$');
h(4)=semilogx(xlim,xlim.*lRatR500(1)+lRatR500(2),'-b','Linewidth',2,'DisplayName','fit');
h(5)=semilogx(xlim,xlim.*lRatM500(1)+lRatM500(2),'-r','Linewidth',2,'DisplayName','fit');


hold off
xlabelmine('$M_{\mathrm{vir}}\,\mathrm{10^{14}M_\odot}]$')
ylabelmine('Ratio');%            '$M_{200},M_{500}\,[\mathrm{M_\odot}]$')
grid
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEastOutside');

set(gca,'fontsize',12,'box','on')


%% 
hf6=figure;

h=[];
h(1)=plot(ll,ratR200,'ob','DisplayName','$R_{200}/R_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);
hold on
h(2)=plot(ll,ratM200,'dr','DisplayName','$M_{200}/M_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);


xlim([0.5 16.5])
ylim([1.02 1.32])

% 
h(3)=plot(xlim,a200.*ones(size(xlim)),'--k','Linewidth',2,'DisplayName','$\sqrt{200/\Delta_{\mathrm{vir}}}$');
h(4)=plot(xlim,xlim.*lRatR200(1)+lRatR200(2),'-b','Linewidth',2,'DisplayName','fit');
h(5)=plot(xlim,xlim.*lRatM200(1)+lRatM200(2),'-r','Linewidth',2,'DisplayName','fit');

hold off
xlabelmine('$R_{\mathrm{vir}}\,[\mathrm{Mpc}]$')
ylabelmine('Ratio');%            '$M_{200},M_{500}\,[\mathrm{M_\odot}]$')
grid
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEastOutside');

set(gca,'fontsize',12,'box','on')

%%
hf7=figure;

h=[];

h(1)=plot(ll,ratR500,'ob','DisplayName','$R_{500}/R_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);
hold on
h(2)=plot(ll,ratM500,'dr','DisplayName','$M_{500}/M_{\mathrm{vir}}$','MarkerSize',10,'Linewidth',2);


xlim([0.5 16.5])
ylim([0.77 0.97])

% 
h(3)=plot(xlim,a500.*ones(size(xlim)),'--k','Linewidth',2,'DisplayName','$\sqrt{500/\Delta_{\mathrm{vir}}}$');
h(4)=plot(xlim,xlim.*lRatR500(1)+lRatR500(2),'-b','Linewidth',2,'DisplayName','fit');
h(5)=plot(xlim,xlim.*lRatM500(1)+lRatM500(2),'-r','Linewidth',2,'DisplayName','fit');


hold off
xlabelmine('$R_{\mathrm{vir}}\,[\mathrm{Mpc}]$')
ylabelmine('Ratio');%            '$M_{200},M_{500}\,[\mathrm{M_\odot}]$')
grid
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEastOutside');

set(gca,'fontsize',12,'box','on')