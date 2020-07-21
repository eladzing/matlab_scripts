
%%  plot evolution comparison

setEnv_RPS


cmap=varycolor(11);% brewermap(11,'Paired');

%depletion ; tim

tdepl=(galResults(11).gasMass(1)./galResults(11).sfr(1))./1e9;


printTag='m11_fgs1';



%% stellar mass 

fgs=1;
figure('position',[1432 421 1000 750],'Color','w')
h=[];

% time 

for i=1:length(galResults)
    
    tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
h(i)=plot(galResults(i).time,galResults(i).stellarMass./galResults(i).stellarMass(1),...
    'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
hold on
end
%plot(galResults(1).time,(1+fgs).*ones(size(galResults(1).time)),':k')

% text(5,1.05,'Max possible growth=1.25')
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',16,'box','on','location','southEast');

%ylim([1 (1+fgs.*1.1)])

xlabelmine('time [Gyr]',16);
ylabelmine('Stellar Mass Growth',16);

set(gca,'Fontsize',14)

if printFlag
    fname=['galEvol_stellarMass_time_' printTag];
    printout_fig(gcf,fname,'v');
end


% rposition

figure('position',[1432 421 1000 750],'Color','w')
h=[];
for i=1:length(galResults)
    
    tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
h(i)=plot(galResults(i).rpos./host.Rvir,galResults(i).stellarMass./galResults(i).stellarMass(1),...
    'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
hold on
end
%plot(galResults(1).rpos,(1+fgs).*ones(size(galResults(1).time)),':k')
text(500,1.05,'Max possible growth=1.25')
grid
xlim([0 1.6])
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',16,'box','on','location','southwest');

%ylim([1 (1+fgs.*1.1)])

xlabelmine('cluster-centric distance $[R_\mathrm{vir}]$',16);
ylabelmine('Stellar Mass Growth',16);
set(gca,'Fontsize',14)

if printFlag
    fname=['galEvol_stellarMass_rpos_' printTag];
    printout_fig(gcf,fname,'v');
end


% % rhoICM 
% 
% figure('position',[1432 421 1000 750],'Color','w')
% h=[];
% for i=1:length(galResults)
%     
%     tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
% h(i)=semilogx(galResults(i).rhoICM,galResults(i).stellarMass./galResults(i).stellarMass(1),...
%     'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
% hold on
% end
% %plot(galResults(1).rpos,(1+fgs).*ones(size(galResults(1).time)),':k')
% text(500,1.05,'Max possible growth=1.25')
% grid
% 
% hl=legend(h);
% set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northEast');
% 
% %ylim([1 (1+fgs.*1.1)])
% 
% xlabelmine('$\rho_{\mathrm{ICM}}$');
% ylabelmine('Stellar Mass Growth');







%% gas mass 

figure('position',[1432 421 1000 750],'Color','w')
h=[];
for i=1:length(galResults)
    
    tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
h(i)=plot(galResults(i).time,galResults(i).gasMass./galResults(i).gasMass(1),...
    'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
hold on
end
%plot(galResults(1).time,(1+fgs).*ones(size(galResults(1).time)),':k')
yl=ylim;
plot(tdepl.*[1 1],yl,'--k')

%text(5,1.05,'Max possible growth=1.25')
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northEast');

%ylim([1 (1+fgs.*1.1)])

xlabelmine('time [Gyr]',16);
ylabelmine('Gas Mass Depletion',16);
set(gca,'Fontsize',14)

if printFlag
    fname=['galEvol_gasMass_time_' printTag];
    printout_fig(gcf,fname,'v');
end



figure('position',[1432 421 1000 750],'Color','w')
h=[];
for i=1:length(galResults)
    
    tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
h(i)=semilogy(galResults(i).rpos./host.Rvir,galResults(i).gasMass./galResults(i).gasMass(1),...
    'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
hold on
end
%plot(galResults(1).rpos,(1+fgs).*ones(size(galResults(1).time)),':k')
%text(500,1.05,'Max possible growth=1.25')
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northwest');
xlim([0 1.6])
%ylim([1 (1+fgs.*1.1)])

xlabelmine('cluster-centric distance $[R_\mathrm{vir}]$',16);
ylabelmine('Gas Mass Depletion',16);
set(gca,'Fontsize',14)

if printFlag
    fname=['galEvol_gasMass_rpos_' printTag];
    printout_fig(gcf,fname,'v');
end


%% sfr 

figure('position',[1432 421 1000 750],'Color','w')
h=[];
for i=1:length(galResults)
    
    tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
h(i)=plot(galResults(i).time,galResults(i).sfr,...
    'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
hold on
end
%plot(galResults(1).time,(1+fgs).*ones(size(galResults(1).time)),':k')
yl=ylim;
plot(tdepl.*[1 1],yl,'--k')

%text(5,1.05,'Max possible growth=1.25')
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northEast');

%ylim([1 (1+fgs.*1.1)])

xlabelmine('time [Gyr]',16);
ylabelmine('SFR $[\mathrm{yr^{-1}}]$',16);
set(gca,'Fontsize',14)

if printFlag
    fname=['galEvol_sfr_time_' printTag];
    printout_fig(gcf,fname,'v');
end


figure('position',[1432 421 1000 750],'Color','w')
h=[];
for i=1:length(galResults)
    
    tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
h(i)=semilogy(galResults(i).rpos./host.Rvir,galResults(i).sfr,...
    'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
hold on
end
%plot(galResults(1).rpos,(1+fgs).*ones(size(galResults(1).time)),':k')
%text(500,1.05,'Max possible growth=1.25')
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northwest');
xlim([0 1.6])
%ylim([1 (1+fgs.*1.1)])

xlabelmine('cluster-centric distance $[R_\mathrm{vir}]$',16);
ylabelmine('SFR $[\mathrm{yr^{-1}}]$',16);
set(gca,'Fontsize',14)

if printFlag
    fname=['galEvol_sfr_rpos_' printTag];
    printout_fig(gcf,fname,'v');
end





%% ssfr 

figure('position',[1432 421 1000 750],'Color','w')
h=[];
for i=1:length(galResults)
    
    tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
h(i)=plot(galResults(i).time,galResults(i).ssfr,...
    'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
hold on
end
%plot(galResults(1).time,(1+fgs).*ones(size(galResults(1).time)),':k')
yl=ylim;
plot(tdepl.*[1 1],yl,'--k')

%text(5,1.05,'Max possible growth=1.25')
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northEast');

%ylim([1 (1+fgs.*1.1)])

xlabelmine('time [Gyr]',16);
ylabelmine('sSFR $[\mathrm{yr^{-1}}]$',16);
set(gca,'Fontsize',14)

if printFlag
    fname=['galEvol_ssfr_time_' printTag];
    printout_fig(gcf,fname,'v');
end


figure('position',[1432 421 1000 750],'Color','w')
h=[];
for i=1:length(galResults)
    
    tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
h(i)=semilogy(galResults(i).rpos./host.Rvir,galResults(i).ssfr,...
    'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
hold on
end
%plot(galResults(1).rpos,(1+fgs).*ones(size(galResults(1).time)),':k')
%text(500,1.05,'Max possible growth=1.25')
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northwest');
xlim([0 1.6])
%ylim([1 (1+fgs.*1.1)])

xlabelmine('cluster-centric distance $[R_\mathrm{vir}]$',16);
ylabelmine('sSFR $[\mathrm{yr^{-1}}]$',16);
set(gca,'Fontsize',14)

if printFlag
    fname=['galEvol_ssfr_rpos_' printTag];
    printout_fig(gcf,fname,'v');
end


% figure
% h=[];
% for i=1:length(galResults)
%     
%     tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
% h(i)=loglog(galResults(i).rhoICM,galResults(i).ssfr,...
%     'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
% hold on
% end
% %plot(galResults(1).rpos,(1+fgs).*ones(size(galResults(1).time)),':k')
% %text(500,1.05,'Max possible growth=1.25')
% grid
% 
% hl=legend(h);
% set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northEast');
% %xlim([0 1.6])
% %ylim([1 (1+fgs.*1.1)])
% 
% xlabelmine('$\rho_{\mathrm{ICM}}$');
% ylabelmine('sSFR $[\mathrm{yr^{-1}}]$');
% 
% 
% figure
% h=[];
% for i=1:length(galResults)
%     
%     tag=sprintf('$V_t/V_c=%2.1f$',orbitBank.fac(i));
% h(i)=semilogx(galResults(i).pram,galResults(i).ssfr,...
%     'linewidth',1.5,'color',cmap(i,:),'DisplayName',tag);
% hold on
% end
% %plot(galResults(1).rpos,(1+fgs).*ones(size(galResults(1).time)),':k')
% %text(500,1.05,'Max possible growth=1.25')
% grid
% 
% hl=legend(h);
% set(hl,'Interpreter','latex','fontsize',16,'box','on','location','northEast');
% %xlim([0 1.6])
% %ylim([1 (1+fgs.*1.1)])
% 
% xlabelmine('$P_{\mathrm{ram}}$');
% ylabelmine('sSFR $[\mathrm{yr^{-1}}]$');
% 
% 
