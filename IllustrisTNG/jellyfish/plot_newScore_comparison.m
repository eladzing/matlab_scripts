ll=0.025:0.05:1.025;
cmap=brewermap(5,'Set1');

myFigure('pos',[33 60 2362  1266]);

t=tiledlayout(2,4);

nexttile;
h=[];
h(1)=histogram(scores.scoreHigh,ll,'facealpha',0.9,'linewidth',0.5,...
    'facecolor',cmap(2,:),'DisplayName','High');
hold on
h(2)=histogram(scores.scoreLow,ll,'facealpha',0.3,'linewidth',0.5,...
    'facecolor',cmap(1,:),'DisplayName','Low');
h(3)=histogram(scores.scoreBase,ll,'DisplayStyle','stairs',...
    'edgecolor','k','linewidth',2,'linestyle',':','DisplayName','Base');
legend(h,'Interpreter','latex','fontsize',12,'location','northeast');
set(gca,'yscale','log')
ylim([500 30000])
%%
nexttile;
h=[];
h(1)=histogram(scores.scoreOverall,ll,'facealpha',0.9,'linewidth',0.5,...
    'facecolor',cmap(2,:),'DisplayName','Overall');
hold on
h(2)=histogram(scores.scoreCombined,ll,'facealpha',0.3,'linewidth',0.5,...
    'facecolor',cmap(1,:),'DisplayName','combined');
h(3)=histogram(scores.scoreBase,ll,'DisplayStyle','stairs',...
    'edgecolor','k','linewidth',2,'linestyle',':','DisplayName','Base');
legend(h,'Interpreter','latex','fontsize',12,'location','northeast');
set(gca,'yscale','log')
ylim([500 30000])
%%


nexttile;
h=[];
h(1)=histogram(scores.noDownVote,ll,'facealpha',0.9,'linewidth',0.5,...
    'facecolor',cmap(2,:),'DisplayName','no Down');
hold on
h(2)=histogram(scores.noUpVote,ll,'facealpha',0.3,'linewidth',0.5,...
    'facecolor',cmap(1,:),'DisplayName','no Up');
h(3)=histogram(scores.scoreBase,ll,'DisplayStyle','stairs',...
    'edgecolor','k','linewidth',2,'linestyle',':','DisplayName','Base');
legend(h,'Interpreter','latex','fontsize',12,'location','northeast');
set(gca,'yscale','log')
ylim([500 30000])
%%
nexttile;
h=[];
h(1)=histogram(scores.noBad,ll,'facealpha',0.9,'linewidth',0.5,...
    'facecolor',cmap(2,:),'DisplayName','no Bad');
hold on
% h(2)=histogram(scores.noUpVote,ll,'facealpha',0.3,'linewidth',0.5,...
%      'facecolor',cmap(1,:),'DisplayName','no Up');
h(2)=histogram(scores.scoreBase,ll,'DisplayStyle','stairs',...
    'edgecolor','k','linewidth',2,'linestyle',':','DisplayName','Base');
 legend(h,'Interpreter','latex','fontsize',12,'location','northeast');
set(gca,'yscale','log')
ylim([500 30000])
%%
nexttile;
h=[];
h(1)=histogram(scores.experts5,ll,'facealpha',0.9,'linewidth',0.5,...
    'facecolor',cmap(2,:),'DisplayName','Experts - 5');
hold on
% h(2)=histogram(scores.noUpVote,ll,'facealpha',0.3,'linewidth',0.5,...
%      'facecolor',cmap(1,:),'DisplayName','no Up');
h(2)=histogram(scores.scoreBase,ll,'DisplayStyle','stairs',...
    'edgecolor','k','linewidth',2,'linestyle',':','DisplayName','Base');
 legend(h,'Interpreter','latex','fontsize',12,'location','northeast');
set(gca,'yscale','log')
ylim([500 30000])
%%
nexttile;
h=[];
h(1)=histogram(scores.scoreCombinedExperts5,ll,'facealpha',0.9,'linewidth',0.5,...
    'facecolor',cmap(2,:),'DisplayName','score+ Experts');
hold on
h(2)=histogram(scores.scoreExperts5,ll,'facealpha',0.3,'linewidth',0.5,...
     'facecolor',cmap(1,:),'DisplayName','comb + Experts');
h(3)=histogram(scores.scoreBase,ll,'DisplayStyle','stairs',...
    'edgecolor','k','linewidth',2,'linestyle',':','DisplayName','Base');
 legend(h,'Interpreter','latex','fontsize',12,'location','northeast');
set(gca,'yscale','log')
ylim([500 30000])
%%
nexttile;
h=[];
h(1)=histogram(scores.scoreAll,ll,'facealpha',0.9,'linewidth',0.5,...
    'facecolor',cmap(2,:),'DisplayName','All');
hold on
% h(2)=histogram(scores.noUpVote,ll,'facealpha',0.3,'linewidth',0.5,...
%      'facecolor',cmap(1,:),'DisplayName','no Up');
h(2)=histogram(scores.scoreBase,ll,'DisplayStyle','stairs',...
    'edgecolor','k','linewidth',2,'linestyle',':','DisplayName','Base');
 legend(h,'Interpreter','latex','fontsize',12,'location','northeast');
set(gca,'yscale','log')
ylim([500 30000])

t.XLabel.String='Score' ;
t.XLabel.FontSize=18;
t.XLabel.Interpreter='latex';
t.YLabel.String='No. of Objects' ;
t.YLabel.FontSize=18;
t.YLabel.Interpreter='latex';
t.Padding='compact';
t.TileSpacing='compact';



%%
[bird, binsize, xl,yl]= histogram2d(scores.scoreBase,scores.scoreAll,ones(size(scores.scoreBase)),'len',21);
bbird=squeeze(bird(:,:,1));
myFigure
imagesc(xl,yl,log10(bbird))
hold on
plot([-0.5 21],[-0.5 21],'--k')
plot([0.775 0.775],[0 1],':k')
plot([0 1],[0.775 0.775],':k')
xlabelmine('Base');
ylabelmine('All');
colormap(brewermap(256,'Reds'))
set(gca,'ydir','normal')
hb=colorbar;
barTitle(hb,'log N')
