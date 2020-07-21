r0=0.1:0.1:50;

%cata=cata10;

mg=zeros(length(cata.Ms),length(r0));
sigas=mg;
sfr=mg;
sfr2=sfr;
%f1=((1-exp(-r0).*(1+r0))./r0.^2).^(1-1.4);
ff1=1e6.^1.4/2.5e-4.*(2.*cata.sigma.*cata.fg.*cata.beta.^2).^(-0.4);
ff2=1e6.^1.4/2.5e-4.*1.4.^2.*(cata.sigma.*cata.fg.*cata.beta.^2).^(-0.4);
for i=1:length(cata.Ms)
    
    
    mg(i,:)=cata.Ms(i).*cata.fg(i).*exp_disk_mass(r0,cata.beta(i));
    sigas(i,:)=mg(i,:)./(pi.*(r0.*cata.rd(i)).^2);
    
    %ms(i,:)=cata.Ms(i).*exp_disk_mass(r0,1);
    sfr(i,:)=sfr_proxy_sigma(mg(i,:),r0.*cata.rd(i));
    sfr2(i,:)=sfr(i,:).*(pi.*(cata.rd(i).*r0).^2);
    
    x=r0.*cata.beta(i);
    f1=((1-exp(-x).*(1+x))./x.^2).^(1-1.4);
    tdep1(i,:)=ff1(i).*f1;
    %tdep1(i,:)=1e6.^1.4/2.5e-4.*(2.*cata.sigma(i).*cata.fg(i).*cata.beta(i).^2).^(-0.4).*f1;  
   
    
    f2=(1-exp(-x).*(1+x))./(1-exp(-1.4.*x).*(1+1.4.*x));
     tdep2(i,:)=ff2(i).*f2;
    %tdep2(i,:)=1e6.^1.4/2.5e-4.*1.4.^2.*(cata.sigma(i).*cata.fg(i).*cata.beta(i).^2).^(-0.4).*f2;  
   
    
end

ff1=1e6.^1.4/2.5e-4.*(2.*cata.sigma.*cata.fg.*cata.beta.^2).^(-0.4);
tdep=sigas./sfr;


%% add travel times 
tv=(1/sqrt(G*deltavir(0)*rho_mean(0).*Ms/Mpc^3*4*pi/3))./(yr*1e9);
for i=1:16
dd(i,:)=[clz0(i).dt2 clz0(i).dt1];
ddd(i,:)=[clz06(i).dt2 clz06(i).dt1];
end

indx=find(r0==5);
h=[];
%maskCat=cata.lambdaMask & 0.7.*cata.qmin>1;
maskCat=cata.lambdaMask & cata.qmin>1;
ll=sum(maskCat);
figure 
[hh, xo]=hist(log10(tdep1(maskCat,indx)),100);
h(1)=stairs(xo,100.*hh./ll,'linewidth',2,'color','k',...
'DisplayName','Catalog');

hold on
h(2)=plot((log10(median(dd(:,2)))+9).*[1 1],ylim,'--b','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0)$');

h(3)=plot((log10(median(ddd(:,2)))+9).*[1 1],ylim,'--r','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0.6)$');

h(4)=plot(log10(tv.*1e9).*[1 1],ylim,'--k','linewidth',2,...
    'DisplayName','$R_{\mathrm{vir}}/V_{\mathrm{vir}})_{z=0}$');
set(gca,'Fontsize',14,'box','on')
hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','NorthEast','box','on')


ylabelmine('$\%$')
xlabelmine('$\log(t_{\mathrm{depletion}})\,[\mathrm{Gyr}]$')


figure 

h(1)=stairs(xo,100.*cumsum(hh)./ll,'linewidth',2,'color',[0 0.7 0],...
'DisplayName','Catalog');

hold on
h(2)=plot((log10(median(dd(:,2)))).*[1 1],ylim,'--b','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0)$');

h(3)=plot((log10(median(ddd(:,2)))).*[1 1],ylim,'--r','linewidth',2,...
    'DisplayName','$\overline{t}_{\mathrm{travel}}(z=0.6)$');

h(4)=plot(log10(tv).*[1 1],ylim,'--k','linewidth',2,...
    'DisplayName','$(R_{\mathrm{vir}}/V_{\mathrm{vir}})_{z=0}$');
set(gca,'Fontsize',14,'box','on')
hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast','box','on')
grid
xlim([-0.3 1.5])
ylabelmine('$\%$')
xlabelmine('$\log(t_{\mathrm{depletion}})\,[\mathrm{Gyr}]$')
