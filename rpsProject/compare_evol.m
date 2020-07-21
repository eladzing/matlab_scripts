%% We carry out a parameter survey of galaxy parameters without changeing the total mass: beta, xi and


% perliminaries
units;
setEnv_RPS;


%% prepare sat parameters
fgs=1;
fbs=0.33;

Ms=1e9;
rd=5;
%Mg=fgs.*Ms;
beta=1;
%Mb=fbs.*Ms;
xi=3;
Mh=1e11;
cs=5;


gal=GALAXY('ms',Ms,'rd',rd,'fgs',fgs,'beta',beta,'fbs',fbs,'xi',xi,'Mh',Mh,'cv',cs);


%% prepare host set parameters

mv0=5e14;
cv0=cvir_Mvir(mv0,0);
fg0=0.15;

% set host

host=NFW('mv',mv0,'cc',cv0,'fg',fg0);


%% prepare orbits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/orbitBank1.mat'])


%% evolve galaxy

for i=1:length(orbitBank.orb)
    fprintf('running %i galaxy of %i \n',i,length(orbitBank.orb));
    galResults(i)=galEvolutionMachine(gal,host,orbitBank.orb(i));
end






% pause
%
% %% figures
% tim=radOrb.t(1:ind).*timeunit;
% timLab='Time [Gyr]';
% % plot mstar & mgas
%
% cmap=brewermap(8,'Set1');
%
% h=[];
% figure
% h(1)=plot(tim,mgas(1:ind)./mgas(1),'color',cmap(1,:),'DisplayName','SF Only');
% hold on
% h(2)=plot(tim,mgas2./mgas2(1),'color',cmap(2,:),'DisplayName','SF \& RPS');
%
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',14)
%
%
% grid
% xlabelmine(timLab);
% ylabelmine('Gas Mass Loss');%$[\mathrm{M_\odot}]$')
% set(gca,'Fontsize',14)
%
% printout_fig(gcf,['sfRps_' orbTag '_compare_mgas'],'v')
%
% % ------------------------------
%
% h=[];
% figure
% h(1)=plot(tim,mstar(1:ind)./mstar(1),'color',cmap(1,:),'DisplayName','SF Only');
% hold on
% h(2)=plot(tim,mstar2./mstar2(1),'color',cmap(2,:),'DisplayName','SF \& RPS');
%
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',14,'location','SouthEast')
%
%
% grid
% xlabelmine(timLab);
% ylabelmine('Stellar Mass Growth'); %$[\mathrm{M_\odot}]$')
% set(gca,'Fontsize',14)
%
% printout_fig(gcf,['sfRps_' orbTag '_compare_mstar'],'v')
%
% % ------------------------------
%
% h=[];
% figure
% h(1)=plot(tim,sum(sfr1(:,1:ind),1)./sum(sfr1(:,1),1),'color',cmap(1,:),'DisplayName','SF Only');
% hold on
% h(2)=plot(tim,sum(sfr2,1)./sum(sfr2(:,1),1),'color',cmap(2,:),'DisplayName','SF \& RPS');
%
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',14)
%
%
% grid
% xlabelmine(timLab);
% ylabelmine('SFR Decay'); %$[\mathrm{M_\odot\,yr^{-1}}]$')
% set(gca,'Fontsize',14)
%
% printout_fig(gcf,['sfRps_' orbTag '_compare_sfr'],'v')
%
% % ------------------------------
%
% h=[];
% figure
% h(1)=plot(tim,sum(sfr1(:,1:ind),1)./mstar(1:ind),'color',cmap(1,:),'DisplayName','SF Only');
% hold on
% h(2)=plot(tim,sum(sfr2,1)./mstar2,'color',cmap(2,:),'DisplayName','SF \& RPS');
% plot(tim,ones(size(tim)).*1e-11,'--k')
%
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',14)
%
%
% grid
% xlabelmine(timLab);
% ylabelmine('sSFR $[\mathrm{yr^{-1}}]$')
% set(gca,'Fontsize',14)
%
% printout_fig(gcf,['sfRps_' orbTag '_compare_ssfr'],'v')
%
% % ------------------------------
%
% db=floor(ind/12);
% ii=1:db:ind;
% ii(end)=ind;
%
% figure
% h=[];
% cnt=0;
% for i=ii
%     cnt=cnt+1;
%     tag=sprintf('%3.2f Gyr',tim(i));
%     h(cnt)=plot(gal.gasMass.rr(2:end),gmass(:,i),'Displayname',tag);
%     hold on
% end
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',14)
%
% grid
% xlabelmine('galactic radius [kpc]');
% ylabelmine('Gas Mass $[\mathrm{M_\odot}]$')
% set(gca,'Fontsize',14)
% titlemine('SF Only')
%
% printout_fig(gcf,['sf_' orbTag '_gasProf'],'v')
%
% figure
% h=[];
% cnt=0;
% for i=ii
%         cnt=cnt+1;
%     tag=sprintf('%3.2f Gyr',tim(i));
%     h(cnt)=plot(gal2.gasMass.rr(2:end),gmass2(:,i),'Displayname',tag);
%     hold on
% end
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',14)
%
% grid
% xlabelmine('galactic radius [kpc]');
% ylabelmine('Gas Mass $[\mathrm{M_\odot}]$')
% set(gca,'Fontsize',14)
% titlemine('SF \& RPS')
%
% printout_fig(gcf,['sfRPS_' orbTag '_gasProf'],'v')
%
%
% % -----------------------
%
%
%     figure
% h=[];
% cnt=0;
% for i=ii
%     cnt=cnt+1;
%     tag=sprintf('%3.2f Gyr',tim(i));
%     h(cnt)=plot(gal.gasMass.rr(2:end),sfr1(:,i),'Displayname',tag);
%     hold on
% end
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',14)
%
% grid
% xlabelmine('galactic radius [kpc]');
% ylabelmine('SFR $[\mathrm{M_\odot\,yr^{-1}}]$')
% set(gca,'Fontsize',14)
% titlemine('SF Only')
%
% printout_fig(gcf,['sf_' orbTag '_sfrProf'],'v')
%
% figure
% h=[];
% cnt=0;
% for i=ii
%         cnt=cnt+1;
%     tag=sprintf('%3.2f Gyr',tim(i));
%     h(cnt)=plot(gal2.gasMass.rr(2:end),sfr2(:,i),'Displayname',tag);
%     hold on
% end
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',14)
%
% grid
% xlabelmine('galactic radius [kpc]');
% ylabelmine('SFR $[\mathrm{M_\odot \,yr^{-1}}]$')
% set(gca,'Fontsize',14)
% titlemine('SF \& RPS')
%
% printout_fig(gcf,['sfRPS_' orbTag '_sfrProf'],'v')
%
%
% % -----------------------
%
%
%
%
%
% % plot gas mass profile
%
%
%
%
%
%
%
%
%
%
%
