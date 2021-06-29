%% try out Mo Mao Whit recipe for disk

units;
numIter=15;
GG=G*km^-2*kpc^-1*Ms;


% assumptions
mstar=1e11; % in solarmass
fg=0.0;  % gas to stellar mass ration in th disk
beta=0.5; % ratio of scale radii for stellar and gas components of disk

lambda=lognrnd(log(0.031),0.57,1) % spin parameter 
md=min(lambda,0.05)  %  0.1; %baryon mass in disk (mstar+mgas)/Mv 
%cc=((3/(4*pi*deltavir(0)*rho_mean(0)))^(1/3)*1e3); %in kpc 

% Halo infor 
Mv=mstar*(1+fg)/md
c=cvir_Mvir(Mv,0)
[rvir mvir tvir vvir]=calculate_virials('mv',Mv);
rvir=rvir*1e3; % change to kpc
A=Mv/(log(1+c)-c/(1+c));
B=c/rvir;

% bulge stuff 
BT=0.1;
Mb=BT*mstar*(1+fg)/(1-BT);
mb=Mb/Mv;
xi=100; %ratio of scale radii for stellar and bulge components of disk

jd=md; %fraction of AM in disk 


%rd_0=sqrt(2)*lambda*jd*md^(-4/3)*(1+fg)^(4/3)*cc*mstar^(1/3)
rd_0=sqrt(2)*lambda*jd/md*(1+fg)*rvir;


r0=0.001:0.001:3;
rr=r0.*rvir;
rds=[]; %zeros(length(numIter)+1,1);
rd=rd_0;
rds(1)=rd;
i=1;
while i<numIter
    
    mdisk=mstar.*(1-exp(-1.*rr./rd).*(1+rr./rd)+fg.*(1-exp(-1.*rr./rd.*beta).*(1+rr./rd.*beta)));
    sigma=mstar/(2*pi*rd^2);
    
    % solve adiabatic contraction 
    ri=solve_MMW_NFW(rr,Mv,rvir,c,mstar,rd,fg,beta,Mb,xi);
    %ri=((1-md).*rr+sqrt((rr.*(md-1)).^2+4.*rr.*rvir./mvir.*mdisk))/2;
    %ri2=((1-md).*rr-sqrt((rr.*(md-1)).^2+4.*rr.*rvir./mvir.*mdisk))/2;
    
    % find circular velocity 
    vcdmSq=GG.*((1-md-mb).*A.*(log(1+B.*ri)-B.*ri./(1+B.*ri)))./rr;
    vcbSq=(GG*Mb).*rr./(rd/xi+rr).^2;
    %vcbSq=GG.*(Mb)./rr;
    %vcdmSq=GG.*((1-md)*mvir/rvir).*ri./rr;
    %vcdmSq2=GG.*((1-md)*mvir/rvir).*ri2./rr;
    
    b1=besseli(0,0.5.*rr./rd).*besselk(0,0.5.*rr./rd)-besseli(1,0.5.*rr./rd).*besselk(1,0.5.*rr./rd);
    bBeta=besseli(0,0.5.*rr./rd.*beta).*besselk(0,0.5.*rr./rd.*beta)...
        -besseli(1,0.5.*rr./rd.*beta).*besselk(1,0.5.*rr./rd.*beta);
    vcdkSq=GG.*pi.*sigma.*rd.*(rr./rd).^2.*(b1+fg.*beta^3.*bBeta);
    
    vc=sqrt(vcdmSq+vcdkSq+vcbSq);
    
    % calculate integral 
    eta=rr./rd;
    integrand=eta.^2.*vc.*(exp(-1.*eta)+fg.*beta.^2.*exp(-1.*eta.*beta));
    
    zeta=trapz(eta,integrand)./vvir;
    %rd=sqrt(2)*lambda*jd*md^(-4/3)*(1+fg)^(4/3)*cc*mstar^(1/3)./zeta(end);
    
    rd=sqrt(2)*lambda*jd/md*(1+fg)*rvir./zeta;
    rds(end+1)=rd;
    
    % check to leave loop 
    
    if abs(diff(rds(end-1:end)))/rd < 5e-4  % 
        break
    end
    i=i+1;
end
figure
plot(rds,'.-');
figure
h=[];
h(1)=plot(rr,sqrt(vcdmSq),'-r','DisplayName','DM','LineWidth',2);
hold on
h(2)=plot(rr,sqrt(vcdkSq),'-b','DisplayName','Disk','LineWidth',2);
h(3)=plot(rr,sqrt(vcbSq),'-g','DisplayName','Bulge','LineWidth',2);
h(4)=plot(rr,vc,'-k','DisplayName','Total','LineWidth',2);
h(5)=plot(rd.*[1 1],ylim,':k','DisplayName','$R_d$','LineWidth',2.5);
h(6)=plot(rvir.*[1 1],ylim,'--k','DisplayName','$R_{\mathrm{vir}}$','LineWidth',2);
h(7)=plot(xlim,vvir.*[1 1],'--g','DisplayName','$V_{\mathrm{vir}}$','LineWidth',2);
%ylim([10 200]);
hl=legend(h);
set(hl','Fontsize',14,'Interpreter','latex');
xlabelmine('$r\,[\mathrm{kpc}]$');
ylabelmine('$V_c\,[\mathrm{km/sec}]$');
grid
set(gca,'box','on','FontSize',12)

figure
h=[];
h(1)=plot(rr,sqrt(vcdmSq),'-r','DisplayName','DM','LineWidth',2);
hold on
h(2)=plot(rr,sqrt(vcdkSq),'-b','DisplayName','Disk','LineWidth',2);
h(3)=plot(rr,sqrt(vcbSq),'-g','DisplayName','Bulge','LineWidth',2);
h(4)=plot(rr,vc,'-k','DisplayName','Total','LineWidth',2);
h(5)=plot(rd.*[1 1],ylim,':k','DisplayName','$R_d$','LineWidth',2.5);
h(6)=plot(rvir.*[1 1],ylim,'--k','DisplayName','$R_{\mathrm{vir}}$','LineWidth',2);
h(7)=plot(xlim,vvir.*[1 1],'--g','DisplayName','$V_{\mathrm{vir}}$','LineWidth',2);
xlim([0 25]);
hl=legend(h);
set(hl','Fontsize',14,'Interpreter','latex');
xlabelmine('$r\,[\mathrm{kpc}]$');
ylabelmine('$V_c\,[\mathrm{km/sec}]$');
grid
set(gca,'box','on','FontSize',12)