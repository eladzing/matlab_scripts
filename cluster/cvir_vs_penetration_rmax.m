%% compare cvir to penetration

%% find cvir  - quick and dirty

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

cv=zeros(size(list));
rm1=zeros(size(list));
rm2=zeros(size(list));

cdm=cv;
adm=cv;
aa=cv;
mmv=cv;
ct=cv;
at=cv;


r0=0.01:0.001:1;
units
fac=G*Ms/Mpc/km^2;
for i=1:length(list)
    
    new_env(list(i))
    
    rv=get_rvir;
    mv=get_mvir;
    mmv(i)=mv;
    
    rp=r0.*rv;
    [mg ms mdm] = read_Mass_Profiles(rp);
    mt0= read_MTOT_Profile(rp);
    mt1=mg+ms+mdm;
   
    %ro=mdm./(4*pi/3.*rp.^3)/(deltavir(0)*rho_mean(0));
    %ii=find(ro<1,1,'first');
    %rvp(i)=rp(ii)./rv;
    % rotot2=smooth(rotot,30);
    
    
%     figure(1)
%     h=[];
%     h(1)=loglog(rp./rv,mt0,'b','linewidth',2,'DisplayName','$M_{tot}$');
%     hold on
%     h(2)=loglog(rp./rv,mt1,'k','linewidth',2,'DisplayName','$\Sigma M$');
%     h(3)=loglog(rp./rv,mdm,'r','linewidth',2,'DisplayName','$M_{DM}$');
%     h(4)=loglog(rp./rv,mg,'g','linewidth',2,'DisplayName','$M_{gas}$');
%     h(5)=loglog(rp./rv,ms,'m','linewidth',2,'DisplayName','$M_{star}$');
%     
%     hl=legend(h);
%     set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);
%     
%     xlabelmine('$r/R_{\mathrm{vir}}$')
%     ylabelmine('$M\,[\mathrm{M_\odot}]$')
    
    %% find cvir by fit to NFW 
    [fitresult1, gof1] = fitNFW(rp./rv, mt0./mv);
    [fitresult2, gof2] = fitNFW(rp./rv, mdm./mv);
    
    ct(i)=fitresult1.c;
    at(i)=fitresult1.a;
    
    cdm(i)=fitresult2.c;
    adm(i)=fitresult2.a;
    aa(i)=mdm(end)/mv;
    gf1=gof1.rsquare;
    gf2=gof2.rsquare;
    
    
    %% find rmax
    vc0=sqrt(fac.*mt0./rp);
    vc1=sqrt(fac.*mt1./rp);
    vcdm=sqrt(fac.*mdm./rp);
    vcg=sqrt(fac.*mg./rp);
    vcs=sqrt(fac.*ms./rp);
    
    [~,i1]=max(vc0);
    [~,i2]=max(vcdm);
    
    
    rm1(i)=rp(i1);
    rm2(i)=rp(i2);
    cv(i)=rv./rp(i2);
   %p=100.*(rp(i1)/rp(i2)-1)
    
    
%     figure(2)
%     h=[];
%     h(1)=plot(rp./rv,vc0,'b','linewidth',2,'DisplayName','$V_{tot}$');
%     hold on
%     h(2)=plot(rp./rv,vc1,'k','linewidth',2,'DisplayName','$\Sigma V$');
%     h(3)=plot(rp./rv,vcdm,'r','linewidth',2,'DisplayName','$V_{DM}$');
%     h(4)=plot(rp./rv,vcg,'g','linewidth',2,'DisplayName','$V_{gas}$');
%     h(5)=plot(rp./rv,vcs,'m','linewidth',2,'DisplayName','$V_{star}$');
%     h(6)=plot(rp(i1)./rv.*[1 1],ylim,'--b','linewidth',2,'DisplayName','$R_{max,1}$');
%     h(7)=plot(rp(i2)./rv.*[1 1],ylim,'--r','linewidth',2,'DisplayName','$R_{max,2}$');
%     
%     hl=legend(h);
%     set(hl,'Interpreter','latex','Location','southEast','Fontsize',14);
%     
%     xlabelmine('$r/R_{\mathrm{vir}}$')
%     ylabelmine('$V_c\,[\mathrm{km/sec}]$')
    
    %     [dlro1,lr1]=derive1(log10(rotot),log10(rp));
    %     %[dlro2,lr2]=derive1(log10(rotot-rog),log10(rp));
    %     [dlro2,lr2]=derive1(log10(rotot2'),log10(rp));
    %     dlro3=smooth(dlro1,50);
    %     figure
    %     %plot(lr1-log10(rv),dlro1,lr2-log10(rv),dlro2)
    %     plot(rp(2:end-1)./rv,dlro1,rp(2:end-1)./rv,dlro2,rp(2:end-1)./rv,dlro3')
    %     grid
    %     legend('reg','smooth','smooth2');
    %     xlabelmine('$\log(r)\,[\mathrm{Mpc}]$')
    %     ylabelmine('$\log(\rho)\,[\mathrm{M_\odot/Mpc^3}]$')
    
    %pause
    close all
end

penetration_depth_a1

cc=cv;
figure
subplot(2,2,1)
plot(cc(urlx),maxpen(urlx),'ro')
hold on
plot(cc(rlx),maxpen(rlx),'bo')
ylabelmine('Maximal Penetration $\%$');


subplot(2,2,2)
plot(cc(urlx),maxpen2(urlx),'ro')
hold on
plot(cc(rlx),maxpen2(rlx),'bo')
ylabelmine('Maximal Penetration 2 $\%$');

subplot(2,2,3)
plot(cc(urlx),avgpen(urlx),'ro')
hold on
plot(cc(rlx),avgpen(rlx),'bo')

ylabelmine('Average Penetration $\%$');
xlabelmine('$c_{\mathrm{vir}}$')


subplot(2,2,4)

plot(cc(urlx),streamno(urlx),'ro')
hold on
plot(cc(rlx),streamno(rlx),'bo')

xlabelmine('$c_{\mathrm{vir}}$')
ylabelmine('Stream No. $\%$');



