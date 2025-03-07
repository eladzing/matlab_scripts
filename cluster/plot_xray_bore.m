%% plot cl6 xray photon flux along bore 


%% read  0.5 0.8 dir - X
cc=brewermap(7,'Set1');
units;
if readFlag
    cl=6;
    new_env(cl);
    
    boxx=1;
    
    
    global zred
    global NCELL
    global hub
    
    cellsize=(boxx./hub/NCELL.*Mpc);
    yzA=xrayProj('yz',1);  % read projection 
    %yzI=sum(yzA.data(:,:,inds),3);
    
    
%     %SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3);
%     rhoMean=deltavir(0).*rho_mean(zred);
%     svir=get_tvir./(rhoMean.^(2/3));
%     
%     ro=RHOG(boxx)./rhoMean;
%     tm=T(boxx)./get_tvir;
%     s=S(boxx)./svir;
%     pre=(ro.*tm);
%     zt=ZIa(boxx)+ZII(boxx);
%     
%     [vx vy vz]=get_velocities(boxx);
end
%% run the bore through  I
bpz=-0.12;
bpx=0;
dir='y';plane='yz';
borelen=[-boxx/2 boxx/2];tickjump=0.1;
smook=1; % smoothing length in cells:L smoothing is over box of 2k+1
%colors = distinguishable_colors(8);
for id1=1:length(bpz)
    for id2=1:length(bpx)
        indz=ceil((bpz(id1)-(-boxx./2))./(boxx./NCELL));
        indx=ceil((bpx(id2)-(-boxx./2))./(boxx./NCELL));
        %ind1=bind1-smook:1:bind1+smook;
        %ind2=bind2-smook:1:bind2+smook;
        
        %% make bores
        xrBore=mk_bore_xrayProj(yzA.data.*cellsize,'dir',dir,'xind',indz,'smooth',smook);
        xb=sum(xrBore,3);
        
    end
end


%% plot

%clear ro tm s Mach pre grads
blen=length(xb);
iax=1:blen;
baxis=(boxx./blen).*(iax-0.5)-(boxx./2);
%xl=[-0.45 0];xb=
xl=[0.08 0.38];


%% plot bores
%cfY=[-0.378 -0.25];
%shkY=-0.33;
 cfY=0.19; 
 shkY=[0.13 0.244 0.287];

hf=figure;
h=[];
h(1)=semilogy(baxis,xb,'-','linewidth',2,'color',cc(1,:),...
    'DisplayName','Photon Flux');
%'DisplayName','$\rho_{\mathrm{gas}}/\rho_{\mathrm{vir}}$');
hold on
h(2)=semilogy(cfY(1).*[1 1],[0.01 1e10],'--k','linewidth',2,...
    'DisplayName','CF');

h(3)=semilogy(shkY(1).*[1 1],[0.01 1e10],':k','linewidth',2,...
    'DisplayName','Shock');
h(3)=semilogy(shkY(2).*[1 1],[0.01 1e10],':k','linewidth',2,...
    'DisplayName','Shock');
h(3)=semilogy(shkY(3).*[1 1],[0.01 1e10],':k','linewidth',2,...
    'DisplayName','Shock');

%h(6)=semilogy(shkY(2).*[1 1],[0.01 10.0],':k','linewidth',2,...
%    'DisplayName','Shock');

xlim(xl);
ylim([1e3 1e5])
%ylim([0.2 6])

xlabelmine('$Y\,[\mathrm{Mpc/h}]$')
ylabelmine('$[\mathrm{Counts\,keV^{-1}\,s^{-1}\,cm^{-2}}]$')


hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')%,'ytick',10.^(-1:1))
grid 'minor'



