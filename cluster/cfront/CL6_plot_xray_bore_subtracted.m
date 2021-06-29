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

%% generate map according to profile

% 
% isoFluxFactor=isothermal_xray_counts(rcube./get_rvir,'fc',0.15,'mv',get_mvir);
% oneCell=read_Apec_spectrum;
% flux=sum(oneCell.spec(:,2));
% 
% isoFlux=isoFluxFactor.*flux;


%% run the bore through  I
bpz=-0.12;
bpx=0;
dir='y';plane='yz';
borelen=[-boxx/2 boxx/2];tickjump=0.1;
smook=0; % smoothing length in cells:L smoothing is over box of 2k+1

for id1=1:length(bpz)
    for id2=1:length(bpx)
        indz=ceil((bpz(id1)-(-boxx./2))./(boxx./NCELL));
        indx=ceil((bpx(id2)-(-boxx./2))./(boxx./NCELL));
        %ind1=bind1-smook:1:bind1+smook;
        %ind2=bind2-smook:1:bind2+smook;
        
        %% make bores
         yzI=trapz(yzA.ebins,yzA.data,3).*cellsize;
         xrBore=mk_bore_xrayProj(yzA.data.*cellsize,'dir',dir,'xind',indz,'smooth',smook);
         xb=trapz(yzA.ebins,xrBore,2); %%./(smook*2+1);
%        xb=mk_bore_xrayProj(isoFlux,'dir',dir,'xind',indz,'smooth',smook);
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
ylim([1e1 5e3])
%ylim([0.2 6])

xlabelmine('$Y\,[\mathrm{Mpc/h}]$')
ylabelmine('$[\mathrm{Counts\,cm^{-2}\,s^{-1}}]$')


hl=legend(h);
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')%,'ytick',10.^(-1:1))
grid 'minor'

% %% plot bore with subtraction 
% subtract background 
%bkgn=((10.^fittedmodel.p2).*(rcube(indz,:)/get_rvir).^fittedmodel.p1)';
bkgn=(10.^fittedmodel.A).*(1+(rcube(indz,:)./get_rvir).^2.*fittedmodel.c.^2).^(0.5-3*fittedmodel.b)';

hf=figure;
h=[];
h(1)=semilogy(baxis,xb./bkgn,'-','linewidth',2,'color',cc(1,:),...
    'DisplayName','Photon Flux');
%'DisplayName','$\rho_{\mathrm{gas}}/\rho_{\mathrm{vir}}$');
hold on
%h(4)=semilogy(baxis,bkgn,'-','linewidth',2,'color',cc(4,:),...
%    'DisplayName','Photon Flux');
h(2)=semilogy(cfY(1).*[1 1],[0.01 10],'--k','linewidth',2,...
    'DisplayName','CF');

h(3)=semilogy(shkY(1).*[1 1],[0.01 10],':k','linewidth',2,...
    'DisplayName','Shock');
h(3)=semilogy(shkY(2).*[1 1],[0.01 10],':k','linewidth',2,...
    'DisplayName','Shock');
h(3)=semilogy(shkY(3).*[1 1],[0.01 10],':k','linewidth',2,...
    'DisplayName','Shock');

%h(6)=semilogy(shkY(2).*[1 1],[0.01 10.0],':k','linewidth',2,...
%    'DisplayName','Shock');

xlim(xl);
ylim([0.6 6])
%ylim([0.2 6])

xlabelmine('$Y\,[\mathrm{Mpc/h}]$')
ylabelmine('Flux/Smoothed Profile')

hl=legend(h(2:end));
set(hl,'Interpreter','latex','fontsize',12,'Location','NorthEast');
set(gca,'fontsize',14,'box','on')%,'ytick',10.^(-1:1))
grid 'minor'
