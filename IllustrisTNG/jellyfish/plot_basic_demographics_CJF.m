%% plot the Basic demographics JF and non-JF galaxies of he CJF project

%% perliminaries
global DEFAULT_MATFILE_DIR
global DEFAULT_PRINTOUT_DIR
colors=brewermap(9,'Set1');
%figPos=[1043         181         840         733];
figPos=[1169         181         714         573];
outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/paper'];
%printFlag=true;

 axFont=18;
 legFont=20;
 labFont=20;
 
%% load data
if setupFlag
    % load object table
    load([DEFAULT_MATFILE_DIR '\cosmic_jellyfish_objectTable.mat']);
    
    % load galaxy properties
    load([DEFAULT_MATFILE_DIR '\jf_galProperties_CJF.mat']);
    
    %% add additional data fields
    massRatio=galProps.galStellarMass./galProps.hostM200c;
    zreds=round(illustris.utils.snap2redshift(galProps.snap),2);
    
    %% Define JF
    threshJF=0.8;
    fprintf('JF Threshold set to %i and above. \n', threshJF);
    
    maskJF=objectTable.scoreWeighted>=0.8;
    maskNJF=objectTable.scoreWeighted<=0.2;
    
    % list snaps and redshifts
    snaps=unique(objectTable.snap);
    redshifts=round(illustris.utils.snap2redshift(snaps),2);
    
    %% adress 3 strange objects
    
    mm=galProps.hostM200c<1e11 & maskJF';
    %ind=find(mm); index of 3 weird objects
    %maskJF_full=maskJF;
    %maskJF(mm)=false;
    
    mask100=objectTable.sim=="TNG100";
    mask50=~mask100;
end

%% plot some general info about the sample
if 1==1
    nObject=height(objectTable);
    fprintf('Total number of objects = %i \n',nObject);
    fprintf('Total number of JF = %i (%4.2f %%) \n',sum(maskJF),sum(maskJF)/nObject*100);
    
    
    fprintf('TNG100 number of objects = %i (%4.2f %%) \n',sum(mask100),sum(mask100)/nObject*100);
    fprintf('TNG50 number of objects = %i (%4.2f %%) \n',sum(mask50),sum(mask50)/nObject*100);
    
    fprintf('TNG100 number of JFs = %i (%4.2f %%) \n',sum(maskJF & mask100),sum(maskJF & mask100)/sum(mask100)*100);
    fprintf('TNG50 number of JFs = %i (%4.2f %%) \n',sum(maskJF & mask50),sum(maskJF & mask50)/sum(mask50)*100);
end

%% plot scores

if 1==1
    bins=-0.025:0.05:1.025;
    
    % plot initial scores 
    myFigure('pos',figPos);
    yl=[0.0 0.35];
    t=tiledlayout(2,1);
    nexttile
    histogram(objectTable.scoreRaw,bins,'normalization','probability','facecolor',colors(1,:))
    ylim(yl)
    %yl=ylim;
    hold on
    plot(0.775.*ones(size(yl)),yl,':k' ,'linewidth',1.8)
    legend("Full Sample",'Interpreter','latex','FontSize',legFont)
    %xlabelmine('Score');
    set(gca,'fontsize',axFont,'TickLabelInterpreter','latex')
    %ylabelmine('fraction of populaiton');
    %titlemine('All');
    
    nexttile
    yl2=[0.001 0.5];
    histogram(objectTable.scoreRaw(mask50),bins,'normalization','probability','facecolor',colors(2,:))
    hold on
    histogram(objectTable.scoreRaw(mask100),bins,'normalization','probability','facecolor',colors(5,:))
    
    %yl=ylim;
    ylim(yl2)
    hold on
    plot(0.775.*ones(size(yl2)),yl2,':k'  ,'linewidth',1.8)
    legend(["TNG50","TNG100"],'Interpreter','latex','FontSize',legFont)
    set(gca,'Yscale','log','fontsize',axFont,'Ytick',[0.01 0.1],'TickLabelInterpreter','latex')
    xlabelmine('Score',20);
    %ylabelmine('fraction of populaiton');
    %titlemine('TNG50');
    
    ylabel(t,'Fraction of inspected satellites','fontsize',labFont,'interpreter','latex')
    t.TileSpacing='tight';
    t.Padding='compact';
    
    if printFlag
        fname='cjf_scoreHist_initial';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
    
    
    % plot weighted scores 
    myFigure('pos',figPos);
    yl=[0.005 0.5];
%     t=tiledlayout(2,1);
%     nexttile
    histogram(objectTable.scoreWeighted,bins,'normalization','probability',...
        'facecolor',colors(1,:),'edgecolor','none')
    hold on 
    
    histogram(objectTable.scoreRaw,bins,'normalization','probability',...
        'facecolor','none','edgecolor','k')
    
    ylim(yl)
    hold on
    plot(0.775.*ones(size(yl)),yl,':k' ,'linewidth',1.8)
    legend(["Adjusted Score","Raw Score"] ,'Interpreter','latex','FontSize',legFont)
    %xlabelmine('Score');
    set(gca,'Yscale','log','fontsize',axFont,'Ytick',[0.01 0.1],'TickLabelInterpreter','latex')
    xlabelmine('Score',labFont);
    ylabelmine('Fraction of inspected satellites',labFont);
    set(gca,'fontsize',axFont,'TickLabelInterpreter','latex')
    %ylabelmine('fraction of populaiton');
    %titlemine('All');
    
%     nexttile
%     yl2=[0.001 0.5];
%     histogram(objectTable.scoreWeighted(mask50),bins,'normalization','probability',...
%         'facecolor',colors(2,:),'edgecolor','none')
%     hold on
%     histogram(objectTable.scoreWeighted(mask100),bins,'normalization','probability',...
%         'facecolor',colors(5,:),'edgecolor','none')
%     
%     histogram(objectTable.scoreRaw(mask50),bins,'normalization','probability',...
%         'edgecolor',colors(2,:),'facecolor','none')
%     histogram(objectTable.scoreRaw(mask100),bins,'normalization','probability',...
%         'edgecolor',colors(5,:),'facecolor','none')
%     
%     
%     
%     %yl=ylim;
%     ylim(yl2)
%     hold on
%     plot(0.775.*ones(size(yl2)),yl2,':k'  ,'linewidth',1.8)
%     legend(["TNG50","TNG100"],'Interpreter','latex','FontSize',legFont)
%     set(gca,'Yscale','log','fontsize',axFont,'Ytick',[0.01 0.1])
%     xlabelmine('Score',20);
%     %ylabelmine('fraction of populaiton');
%     %titlemine('TNG50');
%     
%     ylabel(t,'fraction of populaiton','fontsize',labFont,'interpreter','latex')
%     t.TileSpacing='tight';
%     t.Padding='compact';
%     
    if printFlag
        fname='cjf_scoreHist_weighted';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
    
    
end

%% plot comparison matrix between raw and weighted scores 
scR=objectTable.scoreRaw;
scW=objectTable.scoreWeighted;
ng=length(scR);


myFigure;
edj=0:0.2:1;
edj(end)=1.01;
len=length(edj)-1;
cdata=zeros(len,len);
for i=1:len
    for j=1:len
        cdata(i,j)=sum((scW>=edj(i) & scW<edj(i+1)) & ...
            (scR>=edj(j) & scR<edj(j+1)) );
    end
end
 
 cdata=cdata./ng.*100;
 cdata= round(cdata.*100)/100;
%  lab={'0-0.1' '0.1-0.2' '0.2-0.3' '0.3-0.4' '0.4-0.5' ...
%      '0.5-0.6' '0.6-0.7' '0.7-0.8' '0.8-0.9' '0.9-1'};
 lab={'0-0.2' '0.2-0.4' '0.4-0.6' '0.6-0.8' '0.8-1'};
% imagesc(cdata)
% caxis([0 10])
% set(gca,'Ydir','normal','xtick',1:5,'xticklabel',lab,...
%     'ytick',1:5,'yticklabel',lab)
% colormap(cmap)


ht=heatmap(lab,fliplr(lab),flipud(cdata));
caxis([0.01 100 ])
colorbar('off')
%colormap(brewermap(256,'OrRd'))
colormap(brewermap(512,'Blues'))
ht.FontSize=14;
ht.xlabel('Initital Score')
ht.ylabel('Weighted Score')
ht.ColorScaling='log';
ht.ColorLimits=[-3 2.5];
ht.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
ht.NodeChildren(3).XAxis.Label.FontSize=labFont;
ht.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
ht.NodeChildren(3).YAxis.Label.FontSize=labFont;
%ht.NodeChildren(3).Title.Interpreter = 'latex';

 if printFlag
        fname='cjf_scoreComp_weighted_raw';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
    
%%
myFigure;
jfThresh=0.8;
jfCJF=scW>=jfThresh;
jf50=scR>=jfThresh;



jfdata=[sum(jfCJF & ~jf50) , sum(jfCJF & jf50);...
    sum(~jfCJF & ~jf50) , sum(~jfCJF & jf50)]./ng*100;
jfdata= round(jfdata.*100)/100;
ht=heatmap(["no " "yes "],["yes " "no "],jfdata);
ht.FontSize=14;
caxis([0 20])
colorbar('off')
ht.xlabel('Raw')
ht.ylabel('Weighted')
fprintf("total agreement on what is or isn't JF (score>=%s)= %s %% \n",...
   num2str(jfThresh), num2str(jfdata(1,2)+jfdata(2,1)))


[bird, binsize, xl,yl]= histogram2d(scR,scW,ones(size(scR)),'len',21);
bbird=squeeze(bird(:,:,1));
myFigure;
imagesc(xl,yl,log10(bbird))
hold on
plot([-0.5 21]./20,[-0.5 21]./20,'--k')
plot([15.5 15.5]./20,[-0.5 21]./20,':k')
plot([-0.5 21]./20,[15.5 15.5]./20,':k')
xlabelmine('Raw Score',labFont);
ylabelmine('Adjusted Score',labFont);
colormap(brewermap(256,'OrRd'))
set(gca,'ydir','normal','fontsize',axFont,'TickLabelInterpreter','latex')
hb=colorbar;
barTitle(hb,'log N','fontsize',labFont)
set(hb,'fontsize',axFont,'TickLabelInterpreter','latex')
 if printFlag
        fname='cjf_scoreComp_weighted_raw_cmap';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
%% plot stellar mass and host mass functions
if 1==1
    
    
    bins=8.0:0.1:12.5;
    hf=myFigure('pos',figPos);
    axes1 = axes('Parent',hf);
    hold(axes1,'on');
    hs=[];
    hs(1)=histogram(log10(galProps.galStellarMass(mask50)),bins,'facecolor',colors(2,:),...
        'DisplayName',"TNG50");
    hold on
    hs(2)=histogram(log10(galProps.galStellarMass(mask100)),bins,'facecolor',colors(1,:),...
        'DisplayName',"TNG100");
    hs(3)=histogram(log10(galProps.galStellarMass),bins,'Displaystyle','stairs','edgecolor','k',...
        'linewidth',1.5,'DisplayName',"TNG50+TNG100");
    %     hold on
    %     hs(2)=histogram(log10(galProps.galStellarMass(mask50)),bins,'facecolor',colors(2,:),...
    %         'DisplayName',"TNG50");
    set(gca,'fontsize',axFont,'box','on','TickLabelInterpreter','latex');%,'Yscale','log')
    legend(hs,'Interpreter','latex','FontSize',legFont,'numcolumns',2,'box','off','location','northwest')
    xlabelmine('log Satellite Stellar Mass',labFont);
    ylabelmine('No. of inspectd satellites',labFont);
    hold(axes1,'off');
    % Create inset axes
    axes2 = axes('Parent',hf,...
        'Position',[0.592 0.494 0.282 0.289]);
    hold(axes2,'on');
    histogram(log10(galProps.galStellarMass(mask100)),bins,'facecolor',colors(2,:));
    histogram(log10(galProps.galStellarMass(mask50)),bins,'facecolor',colors(1,:));
    histogram(log10(galProps.galStellarMass),bins,'Displaystyle','stairs','edgecolor','k',...
        'linewidth',1.5)
    xlim(axes2,[11 12.3]);
    ylim(axes2,[0.8 250]);
    set(gca,'fontsize',axFont,'box','on','Yscale','log','ytick',[1 10 100],'TickLabelInterpreter','latex')
    if printFlag
        fname='cjf_stellarMassFunction';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
    
    % sstellar mass function for JF 
        bins=8.0:0.1:12.5;
    hf=myFigure('pos',figPos);
    axes1 = axes('Parent',hf);
    hold(axes1,'on');
    hs(1)=histogram(log10(galProps.galStellarMass(maskJF & mask50)),bins,'facecolor',colors(2,:),...
        'DisplayName',"TNG50");
    hold on
    hs(2)=histogram(log10(galProps.galStellarMass(mask50)),bins,...
        'Displaystyle','stairs','edgecolor',colors(2,:),...
        'DisplayName',"TNG50");
    hs(3)=histogram(log10(galProps.galStellarMass(maskJF & mask100)),bins,'facecolor',colors(1,:),...
        'DisplayName',"TNG100");
    hs(4)=histogram(log10(galProps.galStellarMass(mask100)),bins,...
        'Displaystyle','stairs','edgecolor',colors(1,:),...
        'DisplayName',"TNG100");

    %     hs(3)=histogram(log10(galProps.galStellarMass),bins,'Displaystyle','stairs','edgecolor','k',...
%         'linewidth',1.5,'DisplayName',"All satellites");
%     %     hold on
    %     hs(2)=histogram(log10(galProps.galStellarMass(mask50)),bins,'facecolor',colors(2,:),...
    %         'DisplayName',"TNG50");
    set(gca,'fontsize',axFont,'box','on','Yscale','log','TickLabelInterpreter','latex')
    legend(hs,'Interpreter','latex','FontSize',legFont,'numcolumns',3,'box','off')
    xlabelmine('log Stellar Mass',labFont);
    ylabelmine('No. of inspectd satellites',labFont);
    if printFlag
        fname='cjf_stellarMassFunction_JF';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end

    
    
    
    % halo mass function
    [~,ia,~]=unique(galProps.hostTag);
    hmass=galProps.hostM200c(ia);
    hsim=galProps.sim(ia);
    
    bins=10:0.1:15;
    hf=myFigure('pos',figPos);
    hs(1)=histogram(log10(hmass(hsim=="TNG50")),bins,'facecolor',colors(2,:),...
        'DisplayName',"TNG50");
    hold on
    hs(2)=histogram(log10(hmass(hsim=="TNG100")),bins,'facecolor',colors(1,:),...
        'DisplayName',"TNG100");
     hs(3)=histogram(log10(hmass),bins,'Displaystyle','stairs','edgecolor','k',...
        'linewidth',1.5,'DisplayName',"All satellites");
    set(gca,'fontsize',axFont,'box','on','TickLabelInterpreter','latex');%,'Yscale','log')
    %legend(hs,'Interpreter','latex','FontSize',legFont)
    xlabelmine('log Host Mass',labFont);
    ylabelmine('No. of Hosts',labFont);
    % Create inset axes
    axes2 = axes('Parent',hf,...
        'Position',[0.583596638655462,0.605692844677138,0.282,0.289]); % [0.592 0.494 0.282 0.289]);
                    
    hold(axes2,'on');
    histogram(log10(hmass(hsim=="TNG50")),bins,'facecolor',colors(2,:));
    hold on
    histogram(log10(hmass(hsim=="TNG100")),bins,'facecolor',colors(1,:));
    histogram(log10(hmass),bins,'Displaystyle','stairs','edgecolor','k',...
        'linewidth',1.5);
    xlim(axes2,[13.5 14.5]);
    ylim(axes2,[0 100]);
    set(gca,'fontsize',axFont,'box','on','Yscale','log','TickLabelInterpreter','latex')
    %     titlemine("host mass function");
    if printFlag
        fname='cjf_hostMassFunction';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
    
    % no. of sat's
    bins=10:0.1:15;
    hf=myFigure('pos',figPos);
    hs(1)=histogram(log10(galProps.hostM200c(mask50)),bins,'facecolor',colors(2,:),...
        'DisplayName',"TNG50");
    hold on
    hs(2)=histogram(log10(galProps.hostM200c(mask100)),bins,'facecolor',colors(1,:),...
        'DisplayName',"TNG100");
     hs(3)=histogram(log10(galProps.hostM200c),bins,'Displaystyle','stairs','edgecolor','k',...
        'linewidth',1.5,'DisplayName',"All satellites");
    set(gca,'fontsize',axFont,'box','on','TickLabelInterpreter','latex');%,'Yscale','log')
    %legend(hs,'Interpreter','latex','FontSize',legFont)
    ylim([0 4500])
    xlabelmine('log Host Mass',labFont);
    ylabelmine('No. of inspectd satellites',labFont);
    %titlemine("Number of Sat's found in hosts");
    if printFlag
        fname='cjf_satNumber';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
    
    
    bins=-6:0.125:-1;
    hf=myFigure('pos',figPos);
    hs(1)=histogram(log10(massRatio(mask50)),bins,'facecolor',colors(2,:),...
        'DisplayName',"TNG50");
    hold on
    hs(2)=histogram(log10(massRatio(mask100)),bins,'facecolor',colors(1,:),...
        'DisplayName',"TNG100");
     hs(3)=histogram(log10(massRatio),bins,'Displaystyle','stairs','edgecolor','k',...
        'linewidth',1.5,'DisplayName',"All satellites");
    set(gca,'fontsize',axFont,'box','on','TickLabelInterpreter','latex');%,'Yscale','log')
    %legend(hs,'Interpreter','latex','FontSize',legFont)
    xlabelmine('log Satellite Stellar/Host Mass ratio',labFont);
    ylabelmine('No. of inspectd satellites',labFont);
    %titlemine("stellar=to-host mass ratio");
    % Create inset axes
    %     axes2 = axes('Parent',hf,...
    %     'Position',[0.2 0.6 0.23 0.23]);
    %     hold(axes2,'on');
    %     histogram(log10(massRatio(mask50)),bins,'facecolor',colors(2,:));
    %     hold on
    %     histogram(log10(massRatio(mask100)),bins,'facecolor',colors(1,:));
    %     xlim(axes2,[-1.6 -1]);
    %     ylim(axes2,[0 500]);
    %     set(gca,'fontsize',axFont,'box','on','Yscale','log')
    if printFlag
        fname='cjf_massRationHist';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
    %% snap distribution
    
    %     hf=myFigure('pos',figPos);
    %     hs(1)=histogram(galProps.snap(mask100),30.5:1:99.5,'facecolor',colors(2,:),...
    %         'DisplayName',"TNG100");
    %     hold on
    %     hs(2)=histogram(galProps.snap(mask50),30.5:1:99.5,'facecolor',colors(1,:),...
    %         'DisplayName',"TNG50");
    %     yl=ylim;
    %     ind=find(hs(1).Values>0);
    %
    %     for i=1:length(ind)
    %
    %         xx=0.5.*sum(hs(1).BinEdges(ind(i)-1:ind(i)));
    %
    %         yy=hs(1).Values(ind(i))+0.06.*diff(yl);
    %         zz=round(illustris.utils.snap2redshift(xx+1))/100;
    %
    %         text(xx,yy,num2str(zz),'fontsize',18,'interpreter','latex')
    %     end
    %
    %     set(gca,'fontsize',axFont,'box','on');
    %     legend(hs,'Interpreter','latex','FontSize',18,'Location','northwest')
    %     xlabelmine('Snapshot',labFont);
    %     ylabelmine('Number of Satellites',labFont);
    %     %titlemine("snapshot/redshift distribution");
    
    
    % another version
    snp100=unique(galProps.snap(mask100));
    z100=round(illustris.utils.snap2redshift(snp100),2);
    snp50=unique(galProps.snap(mask50));
    z50=round(illustris.utils.snap2redshift(snp50),2);
    hs100=histcounts(galProps.snap(mask100),30.5:1:99.5);
    b100=hs100(hs100>0);
    hs50=histcounts(galProps.snap(mask50),30.5:1:99.5);
    b50=hs50(hs50>0);
    
    hf=myFigure('pos',[1563         428         818         659]);
    
    ax1 = axes('position',[0.1300    0.1100    0.7750    0.7713]);
    h(2)=bar(snp100,b100,0.16,'facecolor',colors(2,:),...
        'Displayname','TNG100');
    hold on;
    h(1)=bar(snp50,b50,1,'facealpha',0.6,'facecolor',colors(1,:),...
        'Displayname','TNG50');
    
    
    
    legend(h,'Interpreter','latex','FontSize',18,'Location','northwest')
    set(gca,'fontsize',axFont,'box','on','TickLabelInterpreter','latex');
    xlabelmine('Snapshot',labFont);
    ylabelmine('No. of inspectd satellites',labFont);
    
    ax2 = axes('Position', get(ax1,'Position'), ...
        'XAxisLocation','top', ...
        'Color','none', ...
        'XColor','k');
    ax2.YAxis.Visible = 'off';
    ax2.XTick=ax1.XTick;
    ax2.XTickLabel={'2' '1.5' '1' '0.7' '0.5' '0.4' '0.3' '0.2' '0.1' '0'};
    ax2.XLim = ax1.XLim;
    ax2.FontSize=axFont;
    xlabelmine('redshift',labFont);
    linkprop([ax1 ax2], 'XLim');
    if printFlag
        fname='cjf_snapDistribution';
        printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
    end
end
    

%% examine samples by snapshots and simulations

% % snapshot bins
% hostMassBins=10.^(10:15);
% stellarMassBins=10.^(8:1:12.5);
% zredBins=[0:0.25:0.75 1.0 1.5 2.1]; %0:0.2001:2.1;
%
% massRatBins=10.^(-6:0.5:-1);
%
%
%
% % bin Indices arrays
% hmBI=discretize(galProps.hostM200c,hostMassBins);
% smBI=discretize(galProps.galStellarMass,stellarMassBins);
% mrBI=discretize(massRatio,massRatBins);
% zrBI=discretize(redshifts,zredBins);
%
% % tags
% simTags=unique(galProps.sim);
% simTags(3)=join([simTags(1) '\&' simTags(2)]);

%% fraction vs. stellar mass in halo mass bins

hostMassBins=10.^(10:15);
stellarMassBins=10.^(8.0:0.3:12.5);
%stellarMassBins100=10.^(9.5:0.3:12.5);
% hmBI=discretize(galProps.hostM200c,hostMassBins);
% smBI=discretize(galProps.galStellarMass,stellarMassBins);

%jellyfish.utils.plot_demographics(maskJF,smBI,stellarMassBins,hmBI,hostMassBins,...
%    'log','label','stellar mass','legtag','$M_\mathrm{host}=$');

jellyfish.utils.plot_demographics_2sims(maskJF,galProps.galStellarMass,stellarMassBins,galProps.hostM200c,hostMassBins,mask50,...
    'log','label','log stellar mass $[\mathrm{M_\odot}]$','legtag','$M_\mathrm{host}=$','cind',5:-1:1);
if printFlag
    fname='cjf_jfFrac_demograf_mstar_mhostBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end
%% fraction vs. host mass in Stellar mass bins

hostMassBins=10.^(10:15);
stellarMassBins=10.^(8:1:12.5);
% hmBI=discretize(galProps.hostM200c,hostMassBins);
% smBI=discretize(galProps.galStellarMass,stellarMassBins);

% jellyfish.utils.plot_demographics(maskJF,hmBI,hostMassBins,smBI,stellarMassBins,...
%     'log','label','host mass','legtag','$M_\mathrm{\ast}=$');

jellyfish.utils.plot_demographics_2sims(maskJF,galProps.hostM200c,hostMassBins,galProps.galStellarMass,stellarMassBins,mask50,...
    'log','label','log host $M_\mathrm{200,c}\,[\mathrm{M_\odot}]$','legtag','$M_\mathrm{\ast}=$','legLoc',{'northwest','northeast'},'flip');

if printFlag
    fname='cjf_jfFrac_demograf_mhost_mstarBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end
%% fraction vs. mass ratio in Stellar mass bins

massRatBins=10.^(-6:0.5:-2);
%mrBI=discretize(massRatio,massRatBins);
stellarMassBins=10.^(8:1:12.5);
%smBI=discretize(galProps.galStellarMass,stellarMassBins);

% jellyfish.utils.plot_demographics(maskJF,mrBI,massRatBins,smBI,stellarMassBins,...
%     'log','label','mass ratio','legtag','$M_\mathrm{\ast}=$');
jellyfish.utils.plot_demographics_2sims(maskJF,massRatio,massRatBins,galProps.galStellarMass,stellarMassBins,mask50,...
    'log','label','mass ratio','legtag','$M_\mathrm{\ast}=$');

if printFlag
    fname='cjf_jfFrac_demograf_mrat_mstarBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end


stellarMassBins=10.^(8:0.3:12.5);
% smBI=discretize(galProps.galStellarMass,stellarMassBins);
massRatBins=10.^(-6:1:-2);
% mrBI=discretize(massRatio,massRatBins);
% jellyfish.utils.plot_demographics(maskJF,smBI,stellarMassBins,mrBI,massRatBins,...
%     'log','label','stellar mass','legtag','$M_\mathrm{sat}/M_\mathrm{h}$');
jellyfish.utils.plot_demographics_2sims(maskJF,galProps.galStellarMass,stellarMassBins,massRatio,massRatBins,mask50,...
    'log','label','log stellar mass $[\mathrm{M_\odot}]$','legtag','$M_\mathrm{sat}/M_\mathrm{h}$');
if printFlag
    fname='cjf_jfFrac_demograf_mstar_mratBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end

hostMassBins=10.^(11:15);
massRatBins=10.^(-6:1:-2);
jellyfish.utils.plot_demographics_2sims(maskJF,galProps.hostM200c,hostMassBins,massRatio,massRatBins,mask50,...
    'log','label','log host $M_\mathrm{200,c}\,[\mathrm{M_\odot}]$','legtag','$M_\mathrm{sat}/M_\mathrm{h}$');
if printFlag
    fname='cjf_jfFrac_demograf_mhost_mratBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end


%% fraction vs. mass rtion in host mass bins

massRatBins=10.^(-6:0.5:-2);
%mrBI=discretize(massRatio,massRatBins);
hostMassBins=10.^(10:15);
%hmBI=discretize(galProps.hostM200c,hostMassBins);

% jellyfish.utils.plot_demographics(maskJF,mrBI,massRatBins,hmBI,hostMassBins,...
%     'log','label','mass ratio','legtag','$M_\mathrm{host}=$');

jellyfish.utils.plot_demographics_2sims(maskJF,massRatio,massRatBins,galProps.hostM200c,hostMassBins,mask50,...
    'log','label','log mass ratio','legtag','$M_\mathrm{host}=$','cind',5:-1:1);
if printFlag
    fname='cjf_jfFrac_demograf_mrat_mhostBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end
%% fraction vs. redshift in host mass bins

hostMassBins=10.^(10:15);
%hmBI=discretize(galProps.hostM200c,hostMassBins);
%zredBins=[0 0.05:0.1:0.55 0.75 1.25 1.75 2.2];
%zx=[0:0.1:0.5 0.7 1 1.5 2];
zredBins=[0 0.15 0.35 0.55 0.75 1.25 1.75 2.2];
zx=    [0.0  0.15  0.5 0.7 1 1.5 2];

%zrBI=discretize(zreds,zredBins);
% jellyfish.utils.plot_demographics(maskJF,zreds,zredBins,galProps.hostM200c,hostMassBins,'xx',zx,...
%     'label','redshift','legtag','$M_\mathrm{host}=$');

% jellyfish.utils.plot_demographics(maskJF2,zreds,zredBins,galProps.hostM200c,hostMassBins,'xx',zx,...
%     'label','redshift','legtag','$M_\mathrm{host}=$');

jellyfish.utils.plot_demographics_2sims(maskJF,zreds,zredBins,galProps.hostM200c,hostMassBins,mask50,'xx',zx,...
    'label','redshift','legtag','$M_\mathrm{host}=$','cind',5:-1:1);
if printFlag
    fname='cjf_jfFrac_demograf_zred_mhostBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end

% jellyfish.utils.plot_demographics_2sims(maskJF2,zrBI,zredBins,hmBI,hostMassBins,mask50,'xx',zx,...
%     'label','redshift','legtag','$M_\mathrm{host}=$');

% jellyfish.utils.plot_demographics(maskJF,zrBI,zredBins,hmBI,hostMassBins,'xx',zx,...
%     'label','redshift','legtag','$M_\mathrm{host}=$','tng50',galProps.sim=="TNG50");
%
% jellyfish.utils.plot_demographics(maskJF,zrBI,zredBins,hmBI,hostMassBins,'xx',zx,...
%     'label','redshift','legtag','$M_\mathrm{host}=$','tng100',galProps.sim=="TNG100");


%% fraction vs. redshift in stellar mass bins

stellarMassBins=10.^(8:1:12.5);

zredBins=[0 0.15 0.35 0.55 0.75 1.25 1.75 2.2];
zx=[0.0  0.15  0.5 0.7 1 1.5 2];

% zredBins=[0 0.05:0.1:0.55 0.75 1.25 1.75 2.2];
% zx=[0:0.1:0.5 0.7 1 1.5 2];

jellyfish.utils.plot_demographics_2sims(maskJF,zreds,zredBins,galProps.galStellarMass,stellarMassBins,mask50,'xx',zx,...
    'label','redshift','legtag','$M_\mathrm{\ast}=$','flip');
if printFlag
    fname='cjf_jfFrac_demograf_zred_mstarBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end


%% fraction vs. redshift in mass ratio  bins

massRatBins=10.^(-6:1:-2);
zredBins=[0 0.05:0.1:0.55 0.75 1.25 1.75 2.2];
zx=[0:0.1:0.5 0.7 1 1.5 2];

jellyfish.utils.plot_demographics_2sims(maskJF,zreds,zredBins,massRatio,massRatBins,mask50,'xx',zx,...
    'label','redshift','legtag','log mass ratio$=$');
if printFlag
    fname='cjf_jfFrac_demograf_zred_mratBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end




%% fraction vs. stellarmass, hostmass mass ration  in redshift bins

%zredBins=[0 0.25 0.501 0.75 1.25 1.75 2.2];
zredBins=[0 0.501  1.25 2.2];
%zx=[0  1 2];

%zrBI=discretize(zreds,zredBins);
zleg=["$z=0-0.5$" "$z=0.5-1$" "$z=1-2$" ];


% stellar mass
stellarMassBins=10.^(8:0.3:12.5);
%smBI=discretize(galProps.galStellarMass,stellarMassBins);

% jellyfish.utils.plot_demographics(maskJF,smBI,stellarMassBins,zrBI,zredBins,...
%     'log','label','stellar mass','legend',zleg);
jellyfish.utils.plot_demographics_2sims(maskJF,galProps.galStellarMass,stellarMassBins,zreds,zredBins,mask50,...
     'log','label','log stellar mass $[\mathrm{M_\odot}]$','legend',zleg,'flip');
if printFlag
    fname='cjf_jfFrac_demograf_mstar_zredBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end

% host mass
hostMassBins=10.^(10:15);
% hmBI=discretize(galProps.hostM200c,hostMassBins);

% jellyfish.utils.plot_demographics(maskJF,hmBI,hostMassBins,zrBI,zredBins,...
%     'log','label','host mass','legend',zleg,'cind',[1:5 7]);
% jellyfish.utils.plot_demographics(maskJF2,hmBI,hostMassBins,zrBI,zredBins,...
%     'log','label','host mass','legend',zleg,'cind',[1:5 7]);

jellyfish.utils.plot_demographics_2sims(maskJF,galProps.hostM200c,hostMassBins,zreds,zredBins,mask50,...
    'log','label','log host $M_\mathrm{200,c}\,[\mathrm{M_\odot}]$','legend',zleg,'cind',[1:5 7],'legLoc',{'northwest','northeast'},'flip');
if printFlag
    fname='cjf_jfFrac_demograf_mhost_zredBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end

% mass ratio
massRatBins=10.^(-6:0.5:-2);
%mrBI=discretize(massRatio,massRatBins);

jellyfish.utils.plot_demographics_2sims(maskJF,massRatio,massRatBins,zreds,zredBins,mask50,...
    'log','label','log mass ratio','legend',zleg,'cind',[1:5 7],'flip');
% jellyfish.utils.plot_demographics_2sims(maskJF,mrBI,massRatBins,zrBI,zredBins,mask50,...
%     'log','label','log mass ratio','legend',zleg,'cind',[1:5 7]);
if printFlag
    fname='cjf_jfFrac_demograf_mrat_zredBin';
    printout_fig(gcf,fname,'pdf','v','printoutdir',outdir);
end


