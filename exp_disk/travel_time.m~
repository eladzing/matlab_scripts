list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
boxx=8.*ones(size(list));
boxx([13 15 16])=4;

units;
load('C:\Users\eladzing\OneDrive\cluster\matlab\mat_files\shockedge.mat')

for k=2:2
    
    switch k
        case 1
            aa='a1';
        case 2
            aa='a06';
    end
    switch aa
        case 'a1'
            shEd=shockedge_a1;
        case 'a06'
            shEd=shockedge_a06;
        otherwise
            error('wronger pants')
    end
    %boxx=8;
    r0=1:256;
    for i=1:length(list)
        
        new_env(list(i),aa);
        
        global NCELL
        global hub
        mv(i)=get_mvir;
        vv(i)=get_vvir;
        rv(i)=get_rvir;
        
        %% creat vr profile
        vr = Vr_full(boxx(i));
        mask=vr>0;
        vr(mask)=0;
        rog=RHOG(boxx(i));
        rog(mask)=0;
        p1 = MAKE_PROFILE_FROM_CUBE(vr.*rog);
        p2=  MAKE_PROFILE_FROM_CUBE(rog);
        
        vrProf=p1./p2;
        rr=r0.*(0.5.*boxx(i)/hub/NCELL);
        
        %% get shockedge
        sh1=shEd{i,4};
        sh2=shEd{i,5};
        sh=cat(2,sh1,sh2);
        shk=min(sh(sh>get_rvir));
        
        %% idnetify path
        i1=find(rr>rv(i),1,'first');
        i2=find(rr>shk,1,'first');
        i3=find(rr>2.*rv(i),1,'first');
        if isempty(i2)
            i2=NCELL;
        end
        if isempty(i3)
            i3=NCELL;
        end
        
        
        dt1(i)=trapz(rr(i1:i2).*Mpc,1./abs(vrProf(i1:i2).*km))./(yr.*1e9);
        dt2(i)=trapz(rr(i1:i3).*Mpc,1./abs(vrProf(i1:i3).*km))./(yr.*1e9);
        
        
        
        cl(i).vr=vrProf;
        cl(i).rr=rr;
        cl(i).rv=get_rvir;
        cl(i).shk=sh;
        cl(i).ind=[i1 i2 i3];
        cl(i).dt1=dt1(i);
        cl(i).dt2=dt2(i);
        cl(i).mv=get_mvir;
        
    end
    
    switch aa
        case 'a1'
            clz0=cl;
        case 'a06'
            clz06=cl;
        otherwise
            error('wronger pants2')
    end
end

figure;

h=[];
for i=1:16
    
    h(1)=loglog(clz0(i).mv,clz0(i).dt2,'ob',...
        'Markersize',12,'linewidth',1.5,...
        'DisplayName','$2R_{\mathrm{vir}} \to 1 R_{\mathrm{vir}},\,z=0$');
    h(2)=loglog(clz06(i).mv,clz06(i).dt2,'or',...
        'Markersize',12,'linewidth',1.5,...
        'DisplayName','$2R_{\mathrm{vir}} \to 1 R_{\mathrm{vir}},\,z=0.6$');
    
    hold on
    h(3)=loglog(clz0(i).mv,clz0(i).dt1,'+b',...
        'Markersize',12,'linewidth',1.5,...
        'DisplayName','$\mathrm{Shock}\to 1 R_{\mathrm{vir}},\,z=0$');
    h(4)=loglog(clz06(i).mv,clz06(i).dt1,'+r',...
        'Markersize',12,'linewidth',1.5,...
        'DisplayName','$\mathrm{Shock}\to 1 R_{\mathrm{vir}},\,z=0.6$');
    
end

tv0=(1/sqrt(G*deltavir(0)*rho_mean(0).*Ms/Mpc^3*4*pi/3))./(yr*1e9);
tv06=(1/sqrt(G*deltavir(0.6)*rho_mean(0.6).*Ms/Mpc^3*4*pi/3))./(yr*1e9);
h(5)=loglog([1e13 1e16],tv0.*[1 1],'--k','linewidth',2,...
    'DisplayName','$(R_{\mathrm{vir}}/V_{\mathrm{vir}})_{z=0}$');
h(6)=loglog([1e13 1e16],tv06.*[1 1],'-.k','linewidth',2,...
    'DisplayName','$(R_{\mathrm{vir}}/V_{\mathrm{vir}})_{z=0.6}$');
%'$\frac{R_{\mathrm{vir}}}{V_{\mathrm{vir}}}\big|_{z=0}$'

xlim([2e13 5e15])
ylim([0.5 100])

set(gca,'Fontsize',14,'box','on')
hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','NorthEast','box','on')

grid
grid minor

xlabelmine('$M_{\mathrm{vir}}\,[\mathrm{M_\odot}]$')
ylabelmine('$t_{\mathrm{travel}}\,[\mathrm{Gyr}]$')


figure
h=[];
h(1)=plot(cl(1).rr./cl(1).rv,cl(1).vr,'color','b','linewidth',2,...
    'DisplayName','CL101');
hold on
plot(cl(1).rr(cl(1).ind(2))./cl(1).rv.*[1 1],[-1000 1000],'--','color','b','linewidth',2,...
    'DisplayName','CL101');

h(2)=plot(cl(8).rr./cl(2).rv,cl(8).vr,'color','r','linewidth',2,...
    'DisplayName','CL3');
plot(cl(8).rr(cl(8).ind(2))./cl(8).rv.*[1 1],[-1000 1000],'--','color','r','linewidth',2,...
    'DisplayName','CL101');


h(3)=plot(cl(15).rr./cl(15).rv,cl(15).vr,'color',[0 0.7 0],'linewidth',2,...
    'DisplayName','CL14');
plot(cl(15).rr(cl(15).ind(2))./cl(15).rv.*[1 1],[-1000 1000],'--','color',[0 0.7 0],'linewidth',2,...
    'DisplayName','CL101');
xlim([0.9 5])
ylim([-600 30]);

set(gca,'Fontsize',14,'box','on')
hl=legend(h);
set(hl,'Fontsize',14','Interpreter','latex','Location','NorthEast','box','on')

grid
grid minor

xlabelmine('$r/R_{\mathrm{vir}}$')
ylabelmine('$V_r\,[\mathrm{km/sec}]$')







%
% for i=1:16
%     i1=cl(i).ind(1);
%
%     i2=find(cl(i).rr>2.*cl(i).rr(i1),1,'first');
%      if isempty(i2)
%        i2=NCELL;
%    end
%     dt2(i)=trapz(cl(i).rr(i1:i2).*Mpc,1./abs(cl(i).vr(i1:i2).*km))./(yr.*1e9);
% end


%tt=(rv.*Mpc)./(vv.*km)./(1e9*yr);
