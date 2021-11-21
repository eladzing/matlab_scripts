%% radial profiles of HI and H2 for TNG100 & TNG50
%% setup
if setupFlag
    binEdges=0.1:0.25:2.2;
    
    %% load important things TNG50
    bp=illustris.set_env(50,'nomount');
    snap=99;
    [subs50,fofs50,subsInfo50]=illustris.loadFofSub(snap);
    
    massThresh=10^8.3;
    satMask50=illustris.infrastructure.generateMask('subs',subs50','fofs',fofs50,'mass',massThresh,'snap',snap,'sats');
    %% load HI H2 values by components and set profiles (3d)
    
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    load([DEFAULT_MATFILE_DIR '\hih2Catalog_snp99_' simDisplayName '.mat'])
    hih2Struct50=hih2Struct;
    
    hprofs50=generate_hih2_profiles_3d(hih2Struct50,fofs50,subs50,'mask',satMask50,'mean200','bins',binEdges);
    
    
    clear hih2Struct
    %% load important things TNG100
    
    bp=illustris.set_env(100,'nomount');
    snap=99;
    [subs100,fofs100,subsInfo100]=illustris.loadFofSub(snap);
    
    massThresh=10^9;
    satMask100=illustris.infrastructure.generateMask('subs',subs100','fofs',fofs100,'mass',massThresh,'snap',snap,'sats');
    %% load HI H2 values by components and set profiles (3d)
    
    
    load([DEFAULT_MATFILE_DIR '\hih2Catalog_snp99_' simDisplayName '.mat'])
    hih2Struct100=hih2Struct;
    
    hprofs100=generate_hih2_profiles_3d(hih2Struct100,fofs100,subs100,'mask',satMask100,'mean200','bins',binEdges);
    
end
%% perliminaries


colors=brewermap(8,'Set1');
cind=[2 5 3 1];
htag={'11-12','12-13','13-14','14-15'};
mtag={'8.3-9','9-10','10-11','11-12'};
modelTag=hih2Struct.Gal.Hmodel;

outDir='C:\Users\eladz\Documents\workProjects\IllustrisTNG\printout\obsComp';

figPos=[ 800          42        1062         954];
%863 93 1482 1236]);
%% galaxy count

% in stellar mass bins
hf=figure('color','w','position',figPos);

%titlemine('Gal');
h=[];
yl=[1 800];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).count(:,j),'-o',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).count(:,j),'--x',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','northeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.5.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('No. of Satetllites',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','satCount','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


% average number of satellites per host
hf=figure('color','w','position',figPos);
titlemine('Average number of satellites per host');

h=[];
yl=[0.0005 0.2];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).count(:,j)./hprofs50.hostNums(j),'-o',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).count(:,j)./hprofs100.hostNums(j),'--x',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.6.*diff(xl),yl(1)+0.5.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('No. of Satetllites',16);
    
end



% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[1 1000];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(j).xMed(:,k),hprofs50.byHostStar4(j).count(:,k),'-o',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(j).xMed(:,k),hprofs100.byHostStar4(j).count(:,k),'--x',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','northeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    %yl=ylim;
    %xl=xlim;
    xlim(xl)
    ylim(yl)
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.5.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('No. of Satetllites',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','satCount','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


%% sSFR
% ssfrmed in stellar mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[1e-13 2e-10];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).ssfrMed(:,j),'-o',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).ssfrMed(:,j),'--x',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.5.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('sSFR [1/yr]',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','ssfr','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

% ssfravg in stellar mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[1e-13 2e-10];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).ssfrAvg(:,j),'-o',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).ssfrAvg(:,j),'--x',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.5.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('sSFR [1/yr]',16);
    
end

%  ssfrmed in host mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[1e-13 2e-10];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(j).xMed(:,k),hprofs50.byHostStar4(j).ssfrMed(:,k),'-o',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(j).xMed(:,k),hprofs100.byHostStar4(j).ssfrMed(:,k),'--x',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','northeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    %yl=ylim;
    %xl=xlim;
    xlim(xl)
    ylim(yl)
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.5.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
     ylabelmine('sSFR [1/yr]',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','ssfr','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end
%% normalized HI mass

% in stellar mass bins
hf=figure('color','w','position',figPos);
titlemine('Gal');
h=[];
yl=[0 0.5];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassMedN(1,:,j),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).Gal.hiMassMedN(1,:,j),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','galHI','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[0 0.5];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(j).xMed(:,k),hprofs50.byHostStar4(j).Gal.hiMassMedN(1,:,k),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(j).xMed(:,k),hprofs100.byHostStar4(j).Gal.hiMassMedN(1,:,k),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    xlim(xl)
    ylim(yl)
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','galHI','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

%% normalized HI mass in CGM


% in stellar mass bins
hf=figure('color','w','position',figPos);
%titlemine('Gal');
h=[];
yl=[0 3];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).CGMall.hiMassMedN(1,:,j),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).CGMall.hiMassMedN(1,:,j),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','cgmHI','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(j).xMed(:,k),hprofs50.byHostStar4(j).CGMall.hiMassMedN(1,:,k),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(j).xMed(:,k),hprofs100.byHostStar4(j).CGMall.hiMassMedN(1,:,k),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','cgmHI','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


%% HI deficiency Gal 

% in stellar mass bins
hf=figure('color','w','position',figPos);
titlemine('Gal');
h=[];
yl=[-0.1 0.2];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassMedD(1,:,j),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).Gal.hiMassMedD(1,:,j),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('HI deficiency',16);
    
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','galHIDef','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[-0.1 0.2];

xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(j).xMed(:,k),hprofs50.byHostStar4(j).Gal.hiMassMedD(1,:,k),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(j).xMed(:,k),hprofs100.byHostStar4(j).Gal.hiMassMedD(1,:,k),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    xlim(xl)
    ylim(yl)
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('HI deficiency',16);
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','galHIDef','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


%% HI deficiency CGM 

% in stellar mass bins
hf=figure('color','w','position',figPos);
titlemine('Gal');
h=[];
yl=[-0.6 1.6];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).CGMall.hiMassMedD(1,:,j),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).CGMall.hiMassMedD(1,:,j),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
     ylabelmine('HI deficiency',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','cgmHIDef','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[-0.6 1.6];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(j).xMed(:,k),hprofs50.byHostStar4(j).CGMall.hiMassMedD(1,:,k),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(j).xMed(:,k),hprofs100.byHostStar4(j).CGMall.hiMassMedD(1,:,k),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    xlim(xl)
    ylim(yl)
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
     ylabelmine('HI deficiency',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','cgmHIDef','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


%% HI normalized deficiency Gal 

% in stellar mass bins
hf=figure('color','w','position',figPos);
titlemine('Gal');
h=[];
yl=[-1 1]*10;
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassMedDN(1,:,j),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).Gal.hiMassMedDN(1,:,j),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','galHIDefN','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[-1 1].*0.2;
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(j).xMed(:,k),hprofs50.byHostStar4(j).Gal.hiMassMedD(1,:,k),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(j).xMed(:,k),hprofs100.byHostStar4(j).Gal.hiMassMedD(1,:,k),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    xlim(xl)
    ylim(yl)
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','galHIDefN','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


%% HI normalized deficiency CGM 

% in stellar mass bins
hf=figure('color','w','position',figPos);
titlemine('Gal');
h=[];
yl=[-1 1]*0.2;
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).CGMall.hiMassMedD(1,:,j),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).CGMall.hiMassMedD(1,:,j),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','cgmHIDefN','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
yl=[-1 1].*0.2;
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(j).xMed(:,k),hprofs50.byHostStar4(j).CGMall.hiMassMedD(1,:,k),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(j).xMed(:,k),hprofs100.byHostStar4(j).CGMall.hiMassMedD(1,:,k),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    
    xlim(xl)
    ylim(yl)
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','cgmHIDefN','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end




%% normalized HI mass vs ssfr

% in stellar mass bins
hf=figure('color','w','position',figPos);
%titlemine('Gal');
h=[];
yl=[0 0.8];
xl=log10([1e-13 5e-9]);
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(k).xsfrMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassSFRMedN(1,:,j),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(k).xsfrMed(:,j),hprofs100.byHostStar4(k).Gal.hiMassSFRMedN(1,:,j),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southwest','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southwest','interpreter','latex');
    end
    
    
    
   xlim(xl);
   ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('sSFR',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end
if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','HIMassBySFR','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(j).xsfrMed(:,k),hprofs50.byHostStar4(j).Gal.hiMassSFRMedN(1,:,k),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(j).xsfrMed(:,k),hprofs100.byHostStar4(j).Gal.hiMassSFRMedN(1,:,k),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southwest','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southwest','interpreter','latex');
    end
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('sSFR',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end


if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','HIMassBySFR','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end




%% normalized HI mass in CGM  vs ssfr
% in stellar mass bins
% in stellar mass bins
hf=figure('color','w','position',figPos);
%titlemine('Gal');
h=[];
yl=([0 1.5]);
xl=log10([1e-13 5e-9]);
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(k).xsfrMed(:,j),hprofs50.byHostStar4(k).CGMall.hiMassSFRMedN(1,:,j),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(k).xsfrMed(:,j),hprofs100.byHostStar4(k).CGMall.hiMassSFRMedN(1,:,j),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southwest','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southwest','interpreter','latex');
    end
    
    
    
   xlim(xl);
   ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('sSFR',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end


if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','cgmHIMassBySFR','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


% in host mass bins
hf=figure('color','w','position',figPos);
h=[];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=semilogy(hprofs50.byHostStar4(j).xsfrMed(:,k),hprofs50.byHostStar4(j).CGMall.hiMassSFRMedN(1,:,k),'-',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        if j==1;hold on;end
        h(end+1)=semilogy(hprofs100.byHostStar4(j).xsfrMed(:,k),hprofs100.byHostStar4(j).CGMall.hiMassSFRMedN(1,:,k),'--',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j});
    end
    
    if k==1
        legend(h(1:2:8),'fontsize',14,'location','southwest','interpreter','latex');
    elseif k==2
        legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southwest','interpreter','latex');
    end
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('sSFR',16);
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$',16);
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','cgmHIMassBySFR','hostMass','stellarMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end


close all

%% scrath


hf=figure('color','w');
h=[];
for k=1:4
    subplot(2,2,k);
    
    
    for j=1:4
        
        h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassMedN(1,:,j),'-',...
            'color',colors(cind(j),:),'DisplayName',htag{j});
        if j==1;hold on;end
        h(end+1)=plot(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).Gal.hiMassMedN(1,:,j),'--',...
            'color',colors(cind(j),:),'DisplayName',htag{j});
        
    end
    
    xlabelmine('$r/R_\mathrm{vir}$');
    ylabelmine('$M_\mathrm{HI}/M_\mathrm{*}$');
    
    
    
end






%% plot by HI host mass


figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs50.byHost.xMed(:,j),hprofs50.byHost.Gal.hiMassMedN(1,:,j),'-',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    h(end+1)=plot(hprofs100.byHost.xMed(:,j),hprofs100.byHost.Gal.hiMassMedN(1,:,j),'--',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    %     h(end+1)=plot(hprofs.byHostStar2(1).xMed(:,j),hprofs.byHostStar2(1).Gal.hiMassMed(1,:,j),'--',...
    %         'color',colors(cind(j),:),'DisplayName',[htag{j} ', low']);
    %     h(end+1)=plot(hprofs.byHostStar2(2).xMed(:,j),hprofs.byHostStar2(2).Gal.hiMassMed(1,:,j),':',...
    %         'color',colors(cind(j),:),'DisplayName',[htag{j} ', high']);
    %
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass / stellar Mass  ');
titlemine('Gal');

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs50.byHost.xMed(:,j),hprofs50.byHost.CGMall.hiMassMedN(1,:,j),'-',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    h(end+1)=plot(hprofs100.byHost.xMed(:,j),hprofs100.byHost.CGMall.hiMassMedN(1,:,j),'--',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    %     h(end+1)=plot(hprofs.byHostStar2(1).xMed(:,j),hprofs.byHostStar2(1).Gal.hiMassMed(1,:,j),'--',...
    %         'color',colors(cind(j),:),'DisplayName',[htag{j} ', low']);
    %     h(end+1)=plot(hprofs.byHostStar2(2).xMed(:,j),hprofs.byHostStar2(2).Gal.hiMassMed(1,:,j),':',...
    %         'color',colors(cind(j),:),'DisplayName',[htag{j} ', high']);
    %
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass  / stellar Mass');
titlemine('CGM');




figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMall.hiMassMed(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMin.hiMassMed(1,:,j),'--',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMout.hiMassMed(1,:,j),':',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    
    
end

legend(h(1:3:12),'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass $[\mathrm{M_\odot}]$');
titlemine('CGM');

% normalized

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.Gal.hiMassMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass/stellar mass ');
titlemine('Gal');

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMall.hiMassMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMin.hiMassMedN(1,:,j),'--',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMout.hiMassMedN(1,:,j),':',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    
    
end

legend(h(1:3:12),'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass / stellar mass');
titlemine('CGM');


%% plot h2 mass



figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.Gal.h2MassMed(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('H2 mass $[\mathrm{M_\odot}]$');
titlemine('H2 Gal');

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMall.h2MassMed(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMin.h2MassMed(1,:,j),'--',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMout.h2MassMed(1,:,j),':',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    
    
end

legend(h(1:3:12),'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('H2 mass $[\mathrm{M_\odot}]$');
titlemine(' H2 CGM');

% normalized

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.Gal.h2MassMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('H2 mass/stellar mass ');
titlemine('H2 Gal');

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMall.h2MassMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMin.h2MassMedN(1,:,j),'--',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMout.h2MassMedN(1,:,j),':',...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    
    
end

legend(h(1:3:12),'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('H2 mass / stellar mass');
titlemine(' H2 CGM');

%% plot offset by host


figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.Gal.hiMassSFRMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    % h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byHost.Gal.hiMassAvgDN(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass offset ');
titlemine('Gal');

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMall.hiMassMedDN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    % h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byHost.CGMall.hiMassAvgDN(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
    
    % h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMin.hiMassMed(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
    %
    % h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMout.hiMassMed(1,:,j),':',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
    
    
end
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass offset ');
titlemine('CGM');


%% plot by host mass vs ssfr

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xsfrMed(:,j),hprofs.byHost.Gal.hiMassSFRMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    % h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byHost.Gal.hiMassAvgDN(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('log sSFR');
ylabelmine('HI mass /stellar mass ');
titlemine('Gal');

h=[];
figure
for j=1:4
    h(end+1)=plot(hprofs.byHost.xsfrMed(:,j),hprofs.byHost.CGMall.hiMassSFRMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    % h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byHost.Gal.hiMassAvgDN(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('log sSFR');
ylabelmine('HI mass /stellar mass ');
titlemine('CGM');


%% plot H2 by host mass vs. SSFR

figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xsfrMed(:,j),hprofs.byHost.Gal.h2MassSFRMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    % h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byHost.Gal.hiMassAvgDN(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('log sSFR');
ylabelmine('H2 mass /stellar mass ');
titlemine(' H2 Gal');

h=[];
figure
for j=1:4
    h(end+1)=plot(hprofs.byHost.xsfrMed(:,j),hprofs.byHost.CGMall.h2MassSFRMedN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    % h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byHost.Gal.hiMassAvgDN(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('log sSFR');
ylabelmine('H2 mass /stellar mass ');
titlemine('H2 CGM');


%% compare to mock observations Gal

% prepare deficiency data
ms=(galHI_star(1).xMedian);
hm=galHI_star(1).yMedian;
ms=ms(~isnan(ms));
hm=hm(~isnan(ms));
msMass=interp1(ms,hm,log10(hih2Struct.galMass));

hiDef=(hih2Struct.Gal.GalHIMass(1,:)./hih2Struct.galMass)./10.^(msMass)-1;

res= generate_profile_from_vantagePoint_Catalog(satStructY,fofs,subs,'mean200','mass',massThresh,'hiDef',hiDef);


%% plot




figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.Gal.hiMassMedDN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    %    h(end+1)=plot(res.rposMed(:,j),res.hiDefMedMed(3,:),'--',...
    %        'color',colors(cind(j),:),'DisplayName',[htag{j} ' Mock']);
    % h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byHost.Gal.hiMassAvgDN(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
end

legend(h,'fontsize',14,'location','northwest','interpreter','latex');
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass offset ');
titlemine('Gal');


%% cgpm


figure
h=[];

for j=1:4
    h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMall.hiMassMedDN(1,:,j),...
        'color',colors(cind(j),:),'DisplayName',htag{j});
    if j==1;hold on;end
    % h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byHost.CGMall.hiMassAvgDN(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
    
    % h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMin.hiMassMed(1,:,j),'--',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
    %
    % h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMout.hiMassMed(1,:,j),':',...
    %     'color',colors(cind(j),:),'DisplayName',htag{j});
    
    
end
xlabelmine('$r/R_\mathrm{vir}$');
ylabelmine('HI mass offset ');
titlemine('CGM');


%% plot by stellar mass

% figure
% h=[];
%
% for j=1:3
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.Gal.hiMassMed(1,:,j),...
%         'color',colors(cind(j),:),'DisplayName',mtag{j});
%     if j==1;hold on;end
% end
%
% legend(h,'fontsize',14,'location','northwest','interpreter','latex');
% xlabelmine('$r/R_\mathrm{vir}$');
% ylabelmine('HI mass $[\mathrm{M_\odot}]$');
% titlemine('Gal');
%
% figure
% h=[];
%
% for j=1:3
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.CGMall.hiMassMed(1,:,j),...
%         'color',colors(cind(j),:),'DisplayName',mtag{j});
%     if j==1;hold on;end
%
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.CGMin.hiMassMed(1,:,j),'--',...
%         'color',colors(cind(j),:),'DisplayName',mtag{j});
%
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.CGMout.hiMassMed(1,:,j),':',...
%         'color',colors(cind(j),:),'DisplayName',mtag{j});
%
%
% end
%
% legend(h(1:3:9),'fontsize',14,'location','northwest','interpreter','latex');
% xlabelmine('$r/R_\mathrm{vir}$');
% ylabelmine('HI mass $[\mathrm{M_\odot}]$');
% titlemine('CGM');
%
% % normalized
%
% figure
% h=[];
%
% for j=1:3
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.Gal.hiMassMedN(1,:,j),...
%         'color',colors(cind(j),:),'DisplayName',mtag{j});
%     if j==1;hold on;end
% end
%
% legend(h,'fontsize',14,'location','northwest','interpreter','latex');
% xlabelmine('$r/R_\mathrm{vir}$');
% ylabelmine('HI mass / stellaer mass');
% titlemine('Gal');
%
% figure
% h=[];
%
% for j=1:3
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.CGMall.hiMassMedN(1,:,j),...
%         'color',colors(cind(j),:),'DisplayName',mtag{j});
%     if j==1;hold on;end
%
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.CGMin.hiMassMedN(1,:,j),'--',...
%         'color',colors(cind(j),:),'DisplayName',mtag{j});
%
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.CGMout.hiMassMedN(1,:,j),':',...
%         'color',colors(cind(j),:),'DisplayName',mtag{j});
%
%
% end
%
% legend(h(1:3:9),'fontsize',14,'location','northwest','interpreter','latex');
% xlabelmine('$r/R_\mathrm{vir}$');
% ylabelmine('HI mass / stellaer mass');
% titlemine('CGM');
%
%
%% plot offset by mass

%
% figure
% h=[];
%
% for j=1:3
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.Gal.hiMassMedN(1,:,j),...
%         'color',colors(cind(j),:),'DisplayName',htag{j});
%     if j==1;hold on;end
%     h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byMass.Gal.hiMassAvgD(1,:,j),'--',...
%         'color',colors(cind(j),:),'DisplayName',htag{j});
% end
%
% legend(h,'fontsize',14,'location','northwest','interpreter','latex');
% xlabelmine('$r/R_\mathrm{vir}$');
% ylabelmine('HI mass offset ');
% titlemine('Gal');
%
% figure
% h=[];
%
% for j=1:3
%     h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byMass.CGMall.hiMassMedDN(1,:,j),...
%         'color',colors(cind(j),:),'DisplayName',htag{j});
%     if j==1;hold on;end
%     h(end+1)=plot(hprofs.byHost.xAvg(:,j),hprofs.byMass.CGMall.hiMassAvgDN(1,:,j),'--',...
%         'color',colors(cind(j),:),'DisplayName',htag{j});
%
%     % h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMin.hiMassMed(1,:,j),'--',...
%     %     'color',colors(cind(j),:),'DisplayName',htag{j});
%     %
%     % h(end+1)=plot(hprofs.byHost.xMed(:,j),hprofs.byHost.CGMout.hiMassMed(1,:,j),':',...
%     %     'color',colors(cind(j),:),'DisplayName',htag{j});
%
%
% end
% xlabelmine('$r/R_\mathrm{vir}$');
% ylabelmine('HI mass offset ');
% titlemine('CGM');