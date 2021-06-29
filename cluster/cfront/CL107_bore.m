%% read  0.5 0.8 dir - X
cc=brewermap(7,'Set1');
if readFlag
    cl=107;
    new_env(cl);
    
    boxx=2;
    
    
    global zred
    global NCELL
    
    
    
    %SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3);
    rhoMean=deltavir(0).*rho_mean(zred);
    svir=get_tvir./(rhoMean.^(2/3));
    
    ro=RHOG(boxx)./rhoMean;
    tm=T(boxx)./get_tvir;
    s=S(boxx)./svir;
    pre=(ro.*tm);
    zt=ZIa(boxx)+ZII(boxx);
    
    [vx vy vz]=get_velocities(boxx);
end
%% run the bore through  I
bp1=0.75; %0.5;  %0.75
bp2=0;
dir='X';plane='xy';
borelen=[-boxx/2 boxx/2];tickjump=0.1;
smook=3; % smoothing length in cells:L smoothing is over box of 2k+1
%colors = distinguishable_colors(8);
for id1=1:length(bp1)
    for id2=1:length(bp2)
        ind1=ceil((bp1(id1)-(-boxx./2))./(boxx./NCELL));
        ind2=ceil((bp2(id2)-(-boxx./2))./(boxx./NCELL));
        %ind1=bind1-smook:1:bind1+smook;
        %ind2=bind2-smook:1:bind2+smook;
        
        %% make bores
        bro=mk_bore(ro,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
        btm=mk_bore(tm,'weight',ro,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
        bs=mk_bore(s,'weight',ro,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
        bpre=mk_bore(pre,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
        bzt=mk_bore(zt,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
        %bzt=mk_bore(zt,'dir','Y','zind',ind1,'xind',ind2);
        
        
        switch lower(plane)
            case{'xy','yx'}
                vper2=vz;
                switch lower(dir)
                    case {'x'}
                        vpar=vx; vper1=vy;
                    case {'y'}
                        vpar=vy; vper1=vx;
                    otherwise
                        error('wrong pants')
                end
                
            case{'xz','zx'}
                vper2=vy;
                switch lower(dir)
                    case {'x'}
                        vpar=vx; vper1=vz;
                    case {'z'}
                        vpar=vz; vper1=vx;
                    otherwise
                        error('wrong pants')
                end
                
            case{'zy','yz'}
                vper2=vx;
                switch lower(dir)
                    case {'z'}
                        vpar=vz; vper1=vy;
                    case {'y'}
                        vpar=vy; vper1=vz;
                    otherwise
                        error('wrong pants')
                end
                
                
        end
        
        vperTot=sqrt(vper1.^2+vper2.^2);
        
        bvpar=mk_bore(vpar./get_vvir,'weight',ro,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
        bvper1=mk_bore(vper1./get_vvir,'weight',ro,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
        bvper2=mk_bore(vper2./get_vvir,'weight',ro,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
        bvperT=mk_bore(vperTot./get_vvir,'weight',ro,'dir','X','yind',ind1,'zind',ind2,'smooth',smook);
    end
end


%% plot

%clear ro tm s Mach pre grads
blen=length(bs);
iax=1:blen;
baxis=(boxx./blen).*(iax-0.5)-(boxx./2);
%xl=[-0.45 0];
xl=[-0.45 0.4];


%% plot bores
%cfY=[-0.378 -0.25];
%shkY=-0.33;
 cfY=[0.09 -0.099];
 shkY=-0.379;

hf=figure;
h=[];
h(1)=semilogy(baxis,bro,'-','linewidth',2,'color',cc(1,:),...
    'DisplayName','$\rho_{\mathrm{gas}}$');
%'DisplayName','$\rho_{\mathrm{gas}}/\rho_{\mathrm{vir}}$');
hold on
h(2)=semilogy(baxis,btm,'-','linewidth',2,'color',cc(2,:),...
    'DisplayName','$T$');
%'DisplayName','$T/T_{\mathrm{vir}}$');

h(3)=semilogy(baxis,bs,'-','color',cc(3,:),'linewidth',2,...
    'DisplayName','$S$');
%'DisplayName','$S/S_{\mathrm{vir}}$');

h(4)=semilogy(baxis,bpre,'-','linewidth',2,'color',cc(4,:),...
    'DisplayName','$P$');
%'DisplayName','$P/P_{\mathrm{vir}}$');

h(5)=semilogy(baxis,bzt./0.1,'-','linewidth',2,'color',cc(5,:),...
    'DisplayName','$Z$');

h(6)=semilogy(cfY(1).*[1 1],[0.01 100.0],'--k','linewidth',2,...
    'DisplayName','CF');

h(6)=semilogy(cfY(2).*[1 1],[0.01 100.0],'--k','linewidth',2,...
    'DisplayName','CF');

h(7)=semilogy(shkY(1).*[1 1],[0.01 100.0],':k','linewidth',2,...
    'DisplayName','Shock');

%h(6)=semilogy(shkY(2).*[1 1],[0.01 10.0],':k','linewidth',2,...
%    'DisplayName','Shock');

xlim(xl);
%ylim([0.1 10])
ylim([0.03 10])


xlabelmine('$X\,[\mathrm{Mpc/h}]$')
ylabelmine('Normalized')

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')
grid



hf=figure;
h=[];
h(1)=plot(baxis,bvpar,'-','linewidth',2,'color',cc(1,:),...
    'DisplayName','$v_{\parallel}$');
hold on
h(2)=plot(baxis,bvper1,'-','linewidth',2,'color',cc(2,:),...
    'DisplayName','$v_{\perp,1}$');

h(3)=plot(baxis,bvper2,'-','color',cc(3,:),'linewidth',2,...
    'DisplayName','$v_{\perp,2}$');

h(4)=plot(baxis,bvperT,'-','color',cc(5,:),'linewidth',2,...
    'DisplayName','$v_{\perp}$');


h(5)=semilogy(cfY(1).*[1 1],[-10 10],'--k','linewidth',2,...
    'DisplayName','CF');

h(5)=semilogy(cfY(2).*[1 1],[-10 10],'--k','linewidth',2,...
    'DisplayName','CF');

h(6)=semilogy(shkY(1).*[1 1],[-10 10],':k','linewidth',2,...
    'DisplayName','Shock');

% h(4)=semilogy(shkY(2).*[1 1],[-10 10],':k','linewidth',2,...
%     'DisplayName','Shock');
grid


xlim(xl);
ylim([-1.5 1.5])
xlabelmine('$X\,[\mathrm{Mpc/h}]$')
ylabelmine('$v/V_{\mathrm{vir}}$')

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')

