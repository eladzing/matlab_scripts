%% plot cl107 xray photon flux along bore 


%% read  0.5 0.8 dir - X
cc=brewermap(7,'Set1');
units;
if readFlag
    cl=107;
    new_env(cl);
    
    boxx=2;
    
    
    global zred
    global NCELL
    global hub
    
    cellsize=(boxx./hub/NCELL.*Mpc);
    xyA=xrayProj('xy',2);  % read projection 
    
end

%%  prepare isothermal sphere output
% get positions
cm=[0 0];
[meshX, meshY] = meshgrid(1:NCELL);

meshX = meshX - (NCELL+1)/2 -cm(1);
meshY = meshY - (NCELL+1)/2 -cm(2);
% Fix Units (to be in Mpc)
meshX = meshX * ((boxx/hub)/NCELL);
meshY = meshY * ((boxx/hub)/NCELL);


rcube=sqrt(meshX.^2+meshY.^2) ; % r cube in Mpc


% 
% %% generate map according to profile
% 
% imaj2=fittedmodel.p2.*(rcube./get_rvir).^fittedmodel.p1;
% 


% 
% isoFluxFactor=isothermal_xray_counts(rcube./get_rvir,'fc',0.15,'mv',get_mvir);
% oneCell=read_Apec_spectrum;
% flux=sum(oneCell.spec(:,2));
% 
% isoFlux=isoFluxFactor.*flux;


%% run the bore through  I
bp1=0.5;  %0.75
bp2=0;
dir='X';plane='xy';
borelen=[-boxx/2 boxx/2];tickjump=0.1;
smook=0; % smoothing length in cells:L smoothing is over box of 2k+1

for id1=1:length(bp1)
    for id2=1:length(bp2)
        ind1=ceil((bp1(id1)-(-boxx./2))./(boxx./NCELL));
        ind2=ceil((bp2(id2)-(-boxx./2))./(boxx./NCELL));
        %indz=ceil((bpz(id1)-(-boxx./2))./(boxx./NCELL));
        %indx=ceil((bpx(id2)-(-boxx./2))./(boxx./NCELL));
        %ind1=bind1-smook:1:bind1+smook;
        %ind2=bind2-smook:1:bind2+smook;
        
        %% make bores
         xyI=trapz(xyA.ebins,xyA.data,3).*cellsize;
         xrBore=mk_bore_xrayProj(xyA.data.*cellsize,'dir',dir,'yind',ind1,'smooth',smook);
         xb=trapz(xyA.ebins,xrBore,2);
%        xb=mk_bore_xrayProj(isoFlux,'dir',dir,'xind',indz,'smooth',smook);
    end
end


%% plot
bkgn=(10.^fittedmodel.A).*(1+(rcube(ind1,:)./get_rvir).^2.*fittedmodel.c.^2).^(0.5-3*fittedmodel.b)';
bkgn1=bkgn;
%clear ro tm s Mach pre grads
blen=length(xb);
iax=1:blen;
baxis=(boxx./blen).*(iax-0.5)-(boxx./2);
xl=[-0.45 0];
%xl=[0.08 0.38];


%% plot bores
cfY=[-0.378 -0.25];
shkY=-0.33;

hf=figure;
h=[];
h(1)=semilogy(baxis,xb./bkgn,'-','linewidth',2,'color',cc(1,:),...
    'DisplayName','Photon Flux');
%'DisplayName','$\rho_{\mathrm{gas}}/\rho_{\mathrm{vir}}$');
hold on
h(2)=semilogy(cfY(1).*[1 1],[0.01 1e3],'--k','linewidth',2,...
    'DisplayName','CF');
h(2)=semilogy(cfY(2).*[1 1],[0.01 1e3],'--k','linewidth',2,...
    'DisplayName','CF');

h(3)=semilogy(shkY(1).*[1 1],[0.01 1e3],':k','linewidth',2,...
    'DisplayName','Shock');

%h(6)=semilogy(shkY(2).*[1 1],[0.01 10.0],':k','linewidth',2,...
%    'DisplayName','Shock');

xlim(xl);
ylim([1e0 1e3])
%ylim([0.2 6])

xlabelmine('$X\,[\mathrm{Mpc/h}]$')
ylabelmine('$[\mathrm{Counts\,cm^{-2}\,s^{-1}}]$')


hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')%,'ytick',10.^(-1:1))
grid 'minor'


%% Top bore 

bp1=0.75;
bp2=0;
dir='X';plane='xy';
borelen=[-boxx/2 boxx/2];tickjump=0.1;
smook=0; % smoothing length in cells:L smoothing is over box of 2k+1

for id1=1:length(bp1)
    for id2=1:length(bp2)
        ind1=ceil((bp1(id1)-(-boxx./2))./(boxx./NCELL));
        ind2=ceil((bp2(id2)-(-boxx./2))./(boxx./NCELL));
        %indz=ceil((bpz(id1)-(-boxx./2))./(boxx./NCELL));
        %indx=ceil((bpx(id2)-(-boxx./2))./(boxx./NCELL));
        %ind1=bind1-smook:1:bind1+smook;
        %ind2=bind2-smook:1:bind2+smook;
        
        %% make bores
         xyI=trapz(xyA.ebins,xyA.data,3).*cellsize;
         xrBore=mk_bore_xrayProj(xyA.data.*cellsize,'dir',dir,'yind',ind1,'smooth',smook);
         xb=trapz(xyA.ebins,xrBore,2);
%        xb=mk_bore_xrayProj(isoFlux,'dir',dir,'xind',indz,'smooth',smook);
    end
end


%% plot
bkgn=(10.^fittedmodel.A).*(1+(rcube(ind1,:)./get_rvir).^2.*fittedmodel.c.^2).^(0.5-3*fittedmodel.b)';
%clear ro tm s Mach pre grads
blen=length(xb);
iax=1:blen;
baxis=(boxx./blen).*(iax-0.5)-(boxx./2);
xl=[-0.45 0.4];
%xl=[0.08 0.38];


%% plot bores
 cfY=[0.09 -0.099];
 shkY=-0.379;

hf=figure;
h=[];
h(1)=semilogy(baxis,xb./bkgn,'-','linewidth',2,'color',cc(1,:),...
    'DisplayName','Photon Flux');
%'DisplayName','$\rho_{\mathrm{gas}}/\rho_{\mathrm{vir}}$');
hold on
h(2)=semilogy(cfY(1).*[1 1],[0.01 1e2],'--k','linewidth',2,...
    'DisplayName','CF');
h(2)=semilogy(cfY(2).*[1 1],[0.01 1e2],'--k','linewidth',2,...
    'DisplayName','CF');

h(3)=semilogy(shkY(1).*[1 1],[0.01 1e2],':k','linewidth',2,...
    'DisplayName','Shock');

%h(6)=semilogy(shkY(2).*[1 1],[0.01 10.0],':k','linewidth',2,...
%    'DisplayName','Shock');

xlim(xl);
ylim([1e0 1e2])
%ylim([0.2 6])

xlabelmine('$x\,[\mathrm{Mpc/h}]$')
ylabelmine('$[\mathrm{Counts\,cm^{-2}\,s^{-1}}]$')


hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')%,'ytick',10.^(-1:1))
grid 'minor'






% %% plot bore with subtraction 

% % prepare isothermal sphere output
% % get positions
% cm=[0 0];
% [meshX, meshY] = meshgrid(1:NCELL);
% 
% meshX = meshX - (NCELL+1)/2 -cm(1);
% meshY = meshY - (NCELL+1)/2 -cm(2);
% % Fix Units (to be in Mpc)
% meshX = meshX * ((boxx/hub)/NCELL);
% meshY = meshY * ((boxx/hub)/NCELL);
% 
% 
% rcube=sqrt(meshX.^2+meshY.^2) ; % r cube in Mpc
% 
% isoFluxFactor=isothermal_xray_counts(rcube./get_rvir,'fc',0.15,'mv',get_mvir);
% oneCell=read_Apec_spectrum;
% flux=sum(oneCell.spec(:,2));
% 
% isoFlux=isoFluxFactor.*flux;
% 
