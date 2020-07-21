%% analyze quenched fraction in clusters

sim='100';
snap=99; %z=0

bp=illustris.set_env(sim);

units; % load general unit structure in cgs.
global illUnits
global DEFAULT_MATFILE_DIR
global simDisplayName


%load([DEFAULT_MATFILE_DIR '/fofs_subs_TNG100_z0.mat']);
subs=illustris.groupcat.loadSubhalos(bp,snap);
fofs=illustris.groupcat.loadHalos(bp,snap);

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
%% set mask
%%
massThresh=10^9; % threshold for *stellar* mass

hostBins=[13 13.5 14 14.5];

galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'sats');


ssfrThresh=1e-11;
%% get data
%%
mass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galMask).*illUnits.massUnit;
ssfr=illustris.utils.calc_ssfr(subs);
ssfr=ssfr(galMask);
hostMass=fofs.Group_M_Crit200(subs.SubhaloGrNr(galMask)+1).*illUnits.massUnit;
%% set position in cluster
%%
satPos=subs.SubhaloPos(:,galMask);
hostPos=fofs.GroupPos(:,subs.SubhaloGrNr(galMask)+1);
hostRv=fofs.Group_R_Crit200(subs.SubhaloGrNr(galMask)+1);

satRad=sqrt(sum((satPos-hostPos).^2,1))./hostRv;
%% some global quenched fractions
%%
qmask=ssfr<=ssfrThresh;
qFracGlobal=sum(qmask)./length(qmask);

%% by host mass:
%%
m=hostMass<10.^13 ;
qm=ssfr(m)<=ssfrThresh;
qFrac13=sum(qm)./length(qm)

m=hostMass>10.^13  & hostMass<10^13.7;
qm=ssfr(m)<=ssfrThresh;
qFrac513=sum(qm)./length(qm)

m=hostMass>10.^13.7 & hostMass<10^14;
qm=ssfr(m)<=ssfrThresh;
qFrac14=sum(qm)./length(qm)

m=hostMass>10.^14  & hostMass<10^14.7;
qm=ssfr(m)<=ssfrThresh;
qFrac514=sum(qm)./length(qm)

m=hostMass>10.^14.7  & hostMass<10^15;
qm=ssfr(m)<=ssfrThresh;
qFrac15=sum(qm)./length(qm)
%% by stellar mass:
%%
m=mass>10.^9  & mass<10^9.7;
qm=ssfr(m)<=ssfrThresh;
qFrac59=sum(qm)./length(qm)

m=mass>10.^9.7 & mass<10^10;
qm=ssfr(m)<=ssfrThresh;
qFrac10=sum(qm)./length(qm)

m=mass>10.^10  & mass<10^10.7;
qm=ssfr(m)<=ssfrThresh;
qFrac510=sum(qm)./length(qm)

m=mass>10.^10.7  & mass<10^11;
qm=ssfr(m)<=ssfrThresh;
qFrac11=sum(qm)./length(qm)

m=mass>10.^11.7  & mass<10^12;
qm=ssfr(m)<=ssfrThresh;
qFrac511=sum(qm)./length(qm)



%% find fraction by distance host mass and stellar mass 

dr=0.25;
rbin=0.25:dr:2.25;
hmBin=13:0.5:15;
smBin=9:0.5:12;

satRadbin=ceil(satRad./dr);
satRadbin(satRadbin>length(rbin))=length(rbin);

qbinT=zeros(length(smBin)-1,length(hmBin),length(rbin));
nbinT=qbinT;
for k=1:length(smBin)-1
    ms=log10(mass)>smBin(k) & log10(mass)<smBin(k+1);
    
    for j=1:length(hmBin)
        
        if j==1
            mh=log10(hostMass)<hmBin(j);
        else
            mh=log10(hostMass)<hmBin(j) & log10(hostMass)>hmBin(j-1);
        end
        
        for i=1:length(rbin)
            m=satRadbin==i & mh & ms;
            qm=ssfr(m)<ssfrThresh;
            qbinT(k,j,i)=sum(qm);
            nbinT(k,j,i)=sum(m);
        end
        
    end
end




%% plot quenched fraction by stellar mass and distance 

nbin=squeeze(sum(nbinT,2));
qbin=squeeze(sum(qbinT,2))./nbin;
qbin(nbin==0)=0;


cmap=brewermap(8,'Set2');
figure('position',[ 1756 452 805 724])
%yyaxis left

subplot(6,1,1:4)

h=[];
for j=1:length(smBin)-1
    
    str=sprintf('$10^{%s}<M_\\mathrm{s}<10^{%s}$',num2str(smBin(j)),num2str(smBin(j+1)));
    
    h(j)=plot(rbin,qbin(j,:),'-+','color',cmap(j,:),...
        'DisplayName',str);
    hold on
end
    
%hl=legend(h);
%set(hl,'Interpreter','latex','location','NorthEastOutside')
    grid
%xlabelmine('$<r/R_{200}$')
ylabelmine('quenched fraction')
titlemine('By steller mass')


%figure
subplot(6,1,5:6)
%yyaxis left
h=[];
for j=1:length(smBin)-1
            
    str=sprintf('$10^{%s}<M_\\mathrm{s}<10^{%s}$',num2str(smBin(j)),num2str(smBin(j+1)));
    h(j)=semilogy(rbin,nbin(j,:),'-+','color',cmap(j,:),...
        'DisplayName',str);
    hold on
end

hl=legend(h);
set(hl,'Interpreter','latex','location','NorthEastOutside','fontsize',10)
    grid
xlabelmine('$<r/R_{200}$')
ylabelmine('$N$ galaxies')
%titlemine('By steller mass')


printout_fig(gcf,['quenchedFraction_stellarMass_' simDisplayName],'v')

%% plot quenched fraction by host mass

nbin=squeeze(sum(nbinT,1));
qbin=squeeze(sum(qbinT,1))./nbin;
qbin(nbin==0)=0;

cmap=brewermap(8,'Set1');

figure('position',[ 1756 452 805 724])
subplot(6,1,1:4)
h=[];
for j=1:length(hmBin)
    if j==1
    str=sprintf('$M_\\mathrm{h}<10^{%s}$',num2str(hmBin(j)));
    else
        str=sprintf('$10^{%s}<M_\\mathrm{h}<10^{%s}$',num2str(hmBin(j-1)),num2str(hmBin(j)));
    end
    h(j)=plot(rbin,qbin(j,:),'-+','color',cmap(j,:),...
        'DisplayName',str);
    hold on
end

hl=legend(h);
set(hl,'Interpreter','latex')
    grid
%xlabelmine('$<r/R_{200}$')
ylabelmine('quenched fraction')
titlemine('By host mass')
subplot(6,1,5:6)

h=[];
for j=1:length(hmBin)
    if j==1
    str=sprintf('$M_\\mathrm{h}<10^{%s}$',num2str(hmBin(j)));
    else
        str=sprintf('$10^{%s}<M_\\mathrm{h}<10^{%s}$',num2str(hmBin(j-1)),num2str(hmBin(j)));
    end
    h(j)=semilogy(rbin,nbin(j,:),'-+','color',cmap(j,:),...
        'DisplayName',str);
    hold on
end

%hl=legend(h);
%set(hl,'Interpreter','latex')
    grid
xlabelmine('$<r/R_{200}$')
ylabelmine('$N$ galaxies')

printout_fig(gcf,['quenchedFraction_hostMass_' simDisplayName],'v')

%% plot by host/stellar mass and no radius 

% nbin=squeeze(sum(nbinT,3));
% qbin=squeeze(sum(qbinT,3))./nbin;
% qbin(nbin==0)=0;
% 
% figure
% surf(smBin(1:end-1),hmBin,qbin)
% 
% figure
% surf(nbin)
% 
% 
% %% plot mass- ssfr
% %%
% figure 
% cmap=brewermap(256,'*Spectral');
% scatter(log10(mass),log10(ssfr),3,log10(hostMass))
% colormap(cmap);
% xlabelmine('$\log M_\star\,[\mathrm{M_\odot}]$');
% ylabelmine('$\log \mathrm{sSFR}\,[\mathrm{yr^{-1}}] $');
% hb=colorbar;
% grid
% barTitle(hb,'$\log\,M_\mathrm{host}\,[\mathrm{M_\odot}]$');
% %% plot pos in host distribution
% %%
% [y,x]=hist(log10(satRad),20);
% bar(x,log10(y));
% ylabelmine('$\log N$')
% xlabelmine('$\log r_\mathrm{sat}/R_\mathrm{200,host}$')
% titlemine('TNG100 satellite position in Host')