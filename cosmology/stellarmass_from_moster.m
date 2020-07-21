function res=stellarmass_from_moster(mhalo,show,randFlag)
%% Based on Moster et al. 2010, we match a stellar mass to halo mass

%mstar=10.^(8:0.1:12);
%
% showFlag=false;
%
% if exist('show','var')
%     showFlag=strcmp(show,'show') | strcmp(show,'plot');
% end
%
%
if ~exist('randFlag','var')
    randFlag=true;
else
    randFlag=strcmp(randFlag,'rand');
end

%% create SHM based on Moster et al. values
%mh=10.^(10:0.01:16);
a=[1 2];

lm1=11.884;
lm1Ps=lm1+a.*0.03;
lm1Ms=lm1-a.*0.023;

m1=10.^lm1;
m1Ps=10.^lm1Ps;
m1Ms=10.^lm1Ms;

mM0=0.0282;
mM0Ps=mM0+a.*0.00061;
mM0Ms=mM0-a.*0.00053;

beta=1.057;
betaPs=beta+a.*0.054;
betaMs=beta-a.*0.046;

gamma=0.556;
gammaPs=gamma+a.*0.01;
gammaMs=gamma-a.*0.004;


ms=mhalo.*2.*mM0.*( (mhalo./m1).^(-beta)+(mhalo./m1).^(gamma)).^-1;

msM1=mhalo.*2.*mM0Ps(1).*( (mhalo./m1Ms(1)).^(-betaMs(1))+(mhalo./m1Ps(1)).^(gammaMs(1))).^-1;
msP1=mhalo.*2.*mM0Ms(1).*( (mhalo./m1Ps(1)).^(-betaPs(1))+(mhalo./m1Ms(1)).^(gammaPs(1))).^-1;

msM2=mhalo.*2.*mM0Ps(2).*( (mhalo./m1Ms(2)).^(-betaMs(2))+(mhalo./m1Ps(2)).^(gammaMs(2))).^-1;
msP2=mhalo.*2.*mM0Ms(2).*( (mhalo./m1Ps(2)).^(-betaPs(2))+(mhalo./m1Ms(2)).^(gammaPs(2))).^-1;

%% draw values randomly
% sig=max( abs(log10(ms)-log10(msM1)) , ...
%     abs(log10(ms)-log10(msP1)));
% res=10.^(normrnd(log10(ms),sig));
if randFlag
sig=max(abs(ms-msM1) , ...
    abs(ms-msP1));
res=normrnd(ms,sig);
else 
    res=ms;
end
%res.msMoster=ms;

showFlag=false;
if exist('show','var')
    showFlag=strcmp(show,'show');
end



% figure
% loglog(mhalo,ms,mhalo,msM1,mhalo,msP1)
% hold on
% loglog(mhalo,res,'.')
%



%
%
%
% %% find mean and range for given stellar mass
%
% mHaloMean=interp1(ms,mh,mstar);
% mHaloP=interp1(msP1,mh,mstar);
% mHaloM=interp1(msM1,mh,mstar);
%
% mean=log10(mHaloMean);
% sig=max(abs(log10(mHaloP)-mean),abs(log10(mHaloM)-mean));
% mHalo= %                      mean+2.*sig.*randn(size(mstar)));
%
% %md=mstar./mHalo;
%
%
if showFlag
    
    figure
    h(1)=loglog(mhalo,ms,'-b','linewidth',2,'DisplayName','Moster et al. 09');
    hold on
    %loglog(msP1,mh,'--b','linewidth',1)
    %loglog(msM1,mh,'--b','linewidth',1)
    h(2)=loglog(mhalo,msP2,'--b','linewidth',1.5,'DisplayName','$2\sigma$');
    loglog(mhalo,msM2,'--b','linewidth',1.5)
    
    h(3)=loglog(mhalo,res,'+r','DisplayName','realization','markerSize',10,'linewidth',2,...
        'DisplayName','realization');
    
    %xlim([min(mstar) max(mstar)])
   % ylim([min(mstar) max(mstar)])
    grid(gca,'Minor')
    
    hl=legend(h);
    set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast')
    set(gca,'Fontsize',14)
    xlabelmine('$M_{\mathrm{Halo}}\,[\mathrm{M_\odot}]$');
    ylabelmine('$M_{\mathrm{stars}} \,[\mathrm{M_\odot}]$');
    
end

%
% figure
% h(1)=loglog(ms,mh,'-b','linewidth',2,'DisplayName','Moster et al. 09');
% hold on
% %loglog(msP1,mh,'--b','linewidth',1)
% %loglog(msM1,mh,'--b','linewidth',1)
% h(2)=loglog(msP2,mh,'--b','linewidth',1.5,'DisplayName','$2\sigma$');
% loglog(msM2,mh,'--b','linewidth',1.5)
%
% h(3)=loglog(mstar,mHalo,'+r','DisplayName','realization','markerSize',10,'linewidth',2);
%
% xlim([min(mstar) max(mstar)])
% grid(gca,'Minor')
%
% hl=legend(h);
% set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest')
% set(gca,'Fontsize',14)
% xlabelmine('$M_{\mathrm{stars}}\,[\mathrm{M_\odot}]$')
% ylabelmine('$M_{\mathrm{halo}}\,[\mathrm{M_\odot}]$')
%
%
%
% %% jus tryin
% figure
%
% h=[];
% ll=cat(2,msP2./mh,fliplr(msM2./mh));
% xl=cat(2,ms,fliplr(ms));
%  h(1)=loglog(ms,(ms./mh),'-b','linewidth',2,'DisplayName','Moster et al. 2010');
% hold on
% h(2)=patch(xl,ll,'b','FaceAlpha',0.1,'DisplayName','$2\sigma$');
% h(1)=loglog(ms,(ms./mh),'-b','linewidth',2,'DisplayName','Moster et al. 2010');
%
%
% %loglog(msP1,mh,'--b','linewidth',1)
% %loglog(msM1,mh,'--b','linewidth',1)
% %h(2)=loglog(msP2,msP2./mh,'--b','linewidth',1.5,'DisplayName','$2\sigma$');
% %loglog(msM2,msM2./mh,'--b','linewidth',1.5)
%
% h(3)=loglog(mstar,mstar./mHalo,'+r','DisplayName','realization','markerSize',10,'linewidth',2);
%
% xlim([min(mstar) max(mstar)])
% grid(gca,'Minor')
% ylim([0.003 0.04])
% hl=legend(h);
% set(hl,'Fontsize',14','Interpreter','latex','Location','SouthWest')
% set(gca,'Fontsize',14)
% xlabelmine('$M_{\mathrm{stars}}\,[\mathrm{M_\odot}]$')
% ylabelmine('$M_{\mathrm{stars}}/M_{\mathrm{halo}}$')
%
%
%
%
% figure
% h(1)=loglog(ms,mh,'-b','linewidth',2,'DisplayName','Moster et al. 2010');
% hold on
% %loglog(msP1,mh,'--b','linewidth',1)
% %loglog(msM1,mh,'--b','linewidth',1)
% h(2)=loglog(msP2,mh,'--b','linewidth',1.5,'DisplayName','$2\sigma$');
% loglog(msM2,mh,'--b','linewidth',1.5)
%
% h(3)=loglog(mstar,mHalo,'+r','DisplayName','realization','markerSize',10,'linewidth',2);
%
% xlim([min(mstar) max(mstar)])
% grid(gca,'Minor')
%
% hl=legend(h);
% set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest')
% set(gca,'Fontsize',14)
% xlabelmine('$M_{\mathrm{stars}}\,[\mathrm{M_\odot}]$')
% ylabelmine('$M_{\mathrm{halo}}\,[\mathrm{M_\odot}]$')
%
%
% end