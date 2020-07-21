%% read

if readFlag
    cl=6;
    new_env(cl);
    
    boxx=1;
    
    
    global zred
    global NCELL
    
    
    
    %SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3);
    rhoMean=deltavir(0).*rho_mean(zred);
    svir=get_tvir./(rhoMean.^(2/3));
    
    ro=RHOG(boxx)./rhoMean;
    tm=T(boxx)./get_tvir;
    s=S(boxx)./svir;
    pre=(ro.*tm);
    
    [vx vy vz]=get_velocities(boxx);
end
%% run the bore through
bp1=-0.12; %[-0.3;-0.25;-0.2;-0.15;-0.1;-0.05;0];
bp2=0;
dir='Y';
borelen=[-boxx/2 boxx/2];tickjump=0.1;
smook=1; % smoothing length in cells:L smoothing is over box of 2k+1
for id1=1:length(bp1)
    for id2=1:length(bp2)
        bind1=ceil((bp1(id1)-(-boxx./2))./(boxx./NCELL));
        bind2=ceil((bp2(id2)-(-boxx./2))./(boxx./NCELL));
        ind1=bind1-smook:1:bind1+smook;
        ind2=bind2-smook:1:bind2+smook;
        
        %% make bores
        bro=mk_bore(ro,'dir','Y','zind',ind1,'xind',ind2);
        btm=mk_bore(tm,'weight',ro,'dir','Y','zind',ind1,'xind',ind2);
        bs=mk_bore(s,'weight',ro,'dir','Y','zind',ind1,'xind',ind2);
        bpre=mk_bore(pre,'dir','Y','zind',ind1,'xind',ind2);
        
        switch dir
            case {'x','X','1','yz','YZ'}
                vpar=vx; vper=sqrt(vy.^2+vz.^2);
            case {'y','Y','2','zx','ZX'}
                vpar=vy; vper=sqrt(vx.^2+vz.^2);
            case {'z','Z','3','xy','XY'}
                vpar=vz; vper=sqrt(vy.^2+vx.^2);
        end
        bvpar=mk_bore(vpar./get_vvir,'weight',ro,'dir','Y','zind',ind1,'xind',ind2);
        bvper=mk_bore(vper./get_vvir,'weight',ro,'dir','Y','zind',ind1,'xind',ind2);
        
    end
end

%% plot

%clear ro tm s Mach pre grads
blen=length(bs);
iax=1:blen;
baxis=(boxx./blen).*(iax-0.5)-(boxx./2);
xl=[0 0.45];


%% plot bores
cfY=0.19;
shkY=[0.244 0.287];


hf=figure;
h=[];
h(1)=semilogy(baxis,bro,'-b','linewidth',2,...
    'DisplayName','$\rho_{\mathrm{gas}}/\rho_{\mathrm{vir}}$');
hold on
h(2)=semilogy(baxis,btm,'-r','linewidth',2,...
    'DisplayName','$T/T_{\mathrm{vir}}$');

h(3)=semilogy(baxis,bs,'-','color',[0 0.75 0],'linewidth',2,...
    'DisplayName','$S/S_{\mathrm{vir}}$');

h(4)=semilogy(baxis,bpre,'-m','linewidth',2,...
    'DisplayName','$P/P_{\mathrm{vir}}$');

h(5)=semilogy(cfY.*[1 1],[0.01 10.0],'--k','linewidth',2,...
    'DisplayName','Cold Front');

h(6)=semilogy(shkY(1).*[1 1],[0.01 10.0],':k','linewidth',2,...
    'DisplayName','Shock');

h(6)=semilogy(shkY(2).*[1 1],[0.01 10.0],':k','linewidth',2,...
    'DisplayName','Shock');

xlim(xl);
ylim([0.25 10.5])

xlabelmine('$Y\,[\mathrm{Mpc/h}]$')
%ylabelmine('')

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')
grid 



hf=figure;
h=[];
h(1)=plot(baxis,bvper,'-b','linewidth',2,...
    'DisplayName','$v_{\parallel}/V_{\mathrm{vir}}$');
hold on
h(2)=plot(baxis,bvpar,'-r','linewidth',2,...
    'DisplayName','$v_{\perp}/V_{\mathrm{vir}}$');
h(3)=semilogy(cfY.*[1 1],[-10 10],'--k','linewidth',2,...
    'DisplayName','Cold Front');

h(4)=semilogy(shkY(1).*[1 1],[-10 10],':k','linewidth',2,...
    'DisplayName','Shock');

h(4)=semilogy(shkY(2).*[1 1],[-10 10],':k','linewidth',2,...
    'DisplayName','Shock');
grid 


xlim(xl);
ylim([-1 1.5])
xlabelmine('$Y\,[\mathrm{Mpc/h}]$')
%ylabelmine('$v/V_{\mathrm{vir}}$')

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')


% %subplot(2,1,1)
% semilogy(baxis,[bro btm bs bpre],'linewidth',3);grid;
% xlim([-0.07 0.44]);ylim([1e-1 1e1]);%;xlim(borelen);ylim([1e-3 1e3]);
% xlabel(sprintf('%s [Mpc/h]',dir),'Fontsize',12);ylabel(sprintf('Profiles'),'Fontsize',12);
% set(gca,'XTick',borelen(1):tickjump:borelen(2),'Fontsize',12)
% %legend('\rho/\rho_{mean}','T/T_{vir}', 'S/S_{vir}','P/(T_{vir}\rho_{mean})','|v|/c_s','|grad(S)|','Location','Southwest');
% legend('\rho/\rho_{mean}','T/T_{vir}', 'S/S_{vir}','P/(T_{vir}\rho_{mean})','Location','Southwest');
% title(sprintf( '%s Bore Profile %s (%0.3g,%0.3g)',clustername,dir,bp1(id1),bp2(id2)),'Fontsize',12);
% if strcmp(pflag,'print')
%     saveas(gcf,sprintf('%s/%s_asl_bore%s_%0.3g-%0.3g_sm%d.png',result_dir,clustername,dir,bp1(id1),bp2(id2),2*smook+1));
% end
%
% figure;
% %subplot(2,1,2)
% plot(baxis,[bvpar bvper],'linewidth',3);grid;
% xlim([-0.07 0.44]);%;ylim([5e-3 30]);xlim(borelen);%;ylim([5e-3 30]);
% xlabel(sprintf('%s [Mpc/h]',dir),'Fontsize',12);ylabel(sprintf('Profiles'),'Fontsize',12);
% set(gca,'XTick',borelen(1):tickjump:borelen(2),'Fontsize',12)
% %subplot(2,4,8)
% legend('v_{par}/V_{vir}','v_{per}/V_{vir}','Location','Southwest');
% %set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0;1 0 1])
% if strcmp(pflag,'print')
%     saveas(gcf,sprintf('%s/%s_asl_vbore%s_%0.3g-%0.3g_sm%d.png',result_dir,clustername,dir,bp1(id1),bp2(id2),2*smook+1));
% end
% %    end
% %end


