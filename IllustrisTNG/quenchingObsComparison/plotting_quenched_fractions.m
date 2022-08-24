%% plot quenched fractions

bp=illustris.set_env(50,'nomount');
massThresh=10^8.3;
[subs,fofs,subsInfo]=illustris.loadFofSub(99);

colors=brewermap(8,'Set1');
global simDisplayName
%% generate the quenched fractions
qfCrit=generate_quenched_fraction(fofs,subs,'massThresh',massThresh,'crit');
qfMean=generate_quenched_fraction(fofs,subs,'massThresh',massThresh,'mean');

%% plot stellar mass - ssfr + queched line

% define stellar mass
galMass= illustris.utils.get_stellar_mass(subs); % stellar mass within 2*rhalf
% define ssfr
ssfrBase=double(illustris.utils.calc_ssfr(subs,'base',1e-15));

galMass=galMass(qfCrit.sampleMask);
ssfr=ssfrBase(qfCrit.sampleMask);

xx=log10(massThresh):0.1:12.5;
ms=interp1(qfCrit.mainSeq.msX,qfCrit.mainSeq.msY,xx,'linear','extrap');


myFigure;
plot(log10(galMass),log10(ssfr),'.','color',colors(2,:));
hold on
plot(qfCrit.mainSeq.msX,qfCrit.mainSeq.msY,'-','linewidth',2);
plot(xx,ms,':','linewidth',1.5,'color',colors(1,:));
plot(xx,ms-1,'--','linewidth',2,'color',colors(1,:));
grid;
set(gca,'fontsize',14);
xlabelmine('log stellar mass');
ylabelmine('log sSFR');

%% plotting perliminaries
cind=[2 5 3 1];
htag={'11-12','12-13','13-14','14-15'};
mtag={'8.3-9','9-10','10-11','11-12'};
figPos=[ 800          42        1062         954];
outDir='C:\Users\eladz\Documents\workProjects\IllustrisTNG\printout\obsComp';

yl=[0 1];
xl=[0 3.7];
%% plot quenched fraction by host mass
for kk=1:2
    
    switch kk
        case 1
            qf=qfCrit;
            tag='crit';
        case 2
            qf=qfMean;
            tag='mean';
    end
    hf=myFigure;
    h=[];
    
    
    for j=1:4
        
        msk=qf.byHost.count(:,j)>0;
        px=[qf.byHost.xMed(msk,j)' fliplr(qf.byHost.xMed(msk,j)')];
        py=[qf.byHost.bsci(msk,1,j)' fliplr(qf.byHost.bsci(msk,2,j)')];
        if j==1;hold on;end
        patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none');
        h(end+1)=plot(qf.byHost.xMed(msk,j),qf.byHost.qfrac(msk,j),'-o',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        
        %         h(end+1)=semilogy(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).ssfrMed(:,j),'--x',...
        %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    legend(h,'fontsize',14,'location','southeast','interpreter','latex');
    %     elseif k==2
    %         legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southeast','interpreter','latex');
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.82.*diff(xl),yl(1)+0.92.*diff(yl),...
        simDisplayName,'interpreter','latex',...
        'fontsize',16);
    %
    %     text(xl(1)+0.05.*diff(xl),yl(1)+0.85.*diff(yl),...
    %         ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
    %         'fontsize',16);
    %
    set(gca,'fontsize',14,'box','on')
    xlabelmine(['$r/R_\mathrm{200,' tag '}$']);
    ylabelmine('quenched fraction');
    
    
    if printFlag
        fname=sprintf('quenchedFrac_byHost_%s_snp99_TNG100',tag);
        printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
    end
    
end





%% plot quenched fraction by stellar and host mass


hf=myFigure('pos',figPos);

h=[];

for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        h(end+1)=plot(qfCrit.byHostStar4(k).xMed(:,j),qfCrit.byHostStar4(k).qfrac(:,j),'-o',...
            'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        if j==1;hold on;end
        %         h(end+1)=semilogy(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).ssfrMed(:,j),'--x',...
        %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h,'fontsize',14,'location','southeast','interpreter','latex');
        %     elseif k==2
        %         legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','southeast','interpreter','latex');
    end
    
    
    
    xlim(xl);
    ylim(yl);
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.85.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('quenched fraction');
    
end

if printFlag
    fname=sprintf('hiProfs_%s_%sSplit_%sBins','ssfr','stellarMass','hostMass');
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

