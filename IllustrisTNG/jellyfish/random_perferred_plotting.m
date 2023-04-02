%% plotting scripts for random vs. perferred viewing angle 

%% load tables 

global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'],'objectTableComp');

global DEFAULT_PRINTOUT_DIR
    outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/paper'];
    
    %%
    axFont=18;
 legFont=20;
 labFont=20;
%% plotting histogram 

hf=myFigure('pos',[1066  231  800  668]);
axes1 = axes('Parent',hf);
hold(axes1,'on');

histogram(objectTableComp.scoreWeightedRand)
hold on 
histogram(objectTableComp.scoreWeightedPref)
plot([0.8 0.8],[0 6000],':k','linewidth',1.5)
legend('Random','Optimized','Interpreter','latex','fontsize',16);
%set(gca,'yscale','log')
box(axes1,'on');
hold(axes1,'off');
set(axes1,'fontsize',axFont)
ylim([0 5e3])
xlabelmine('Score');
ylabelmine('No. of objects',labFont);
%titlemine('weighted Score',labFont);

% Create axes
axes2 = axes('Parent',hf,...
    'Position',[0.434411764705882,0.347985611510792,0.38,0.38]);
hold(axes2,'on');
histogram(objectTableComp.scoreWeightedRand)
hold on 
histogram(objectTableComp.scoreWeightedPref)
plot([0.8 0.8],[1 10000],':k','linewidth',1.8)
ylim([50 5e3])
box(axes2,'on');
hold(axes2,'off');
set(gca,'FontSize',axFont,'YMinorTick','on','YScale','log');


fname='cjf_rand_pref_weight_hist';
if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end



%% difference histogram 
bnn=-1.025:0.05:1.025;
myFigure;
% histogram(objectTableComp.scoreRawPref-objectTableComp.scoreRawRand,-20.5:1:20.5)
% hold on
histogram(objectTableComp.scoreWeightedPref-objectTableComp.scoreWeightedRand,bnn)%,-17.5:0.5:17.5)
grid
set(gca,'yscale','log')
xlim([-1.025 1.025])
%legend('Raw','weighted','Interpreter','latex','fontsize',14);
xlabelmine('Score Optimized - Random');
ylabelmine('No. of objects');



%% 2D histogram 
scR=objectTableComp.scoreWeightedRand;
scP=objectTableComp.scoreWeightedPref;
ng=height(objectTableComp);

[bird, binsize, xl,yl]= histogram2d(scR,scP,ones(size(scR)),'len',21);
bbird=squeeze(bird(:,:,1));

% [bird2, binsize2, xl2,yl2]= histogram2d(scR,scP-scR,ones(size(scR)),'len',[43 21],...
%     'xlim',[0 1],'yylim',[-1 1]);
% bbird2=squeeze(bird2(:,:,1));
% 
% [bird3, binsize3, xl3,yl3]= histogram2d(scP,scP-scR,ones(size(scR)),'len',[42 21],...
%     'xlim',[0 1],'yylim',[-1 1]);
% bbird3=squeeze(bird3(:,:,1));



q1=sum(scP>=0.8 & scR<0.8)/ng*100;
q4=sum(scP<0.8 & scR>=0.8)/ng*100;
q2=sum(scP>=0.8 & scR>=0.8)/ng*100;
q3=sum(scP<0.8 & scR<0.8)/ng*100;

%% plot 2D histogram
hf=myFigure;
imagesc(xl,yl,log10(bbird./sum(sum(bbird)).*100))
hold on
plot([-0.5 21]./20,[-0.5 21]./20,'--k','linewidth',1.5)
plot([15.5 15.5]./20,[-0.5 21]./20,':k','linewidth',1.5)
plot([-0.5 21]./20,[15.5 15.5]./20,':k','linewidth',1.5)
xlabelmine('Random',labFont);
ylabelmine('Optimized',labFont);
%colormap(brewermap(256,'Reds'))
colormap(flipud(magma(256)));
caxis(log10([0.01 20])); 
text(0.05,0.95,[num2str(q1,2) '\%'],'Interpreter','latex',...
    'fontsize',22)
text(0.8,0.95,[num2str(q2,2) '\%'],'Interpreter','latex',...
    'fontsize',22)
text(0.05,0.7,[num2str(q3,2) '\%'],'Interpreter','latex',...
    'fontsize',22)
text(0.8,0.7,[num2str(q4,2) '\%'],'Interpreter','latex',...
    'fontsize',22)

set(gca,'ydir','normal','fontsize',axFont)
hb=colorbar;
barTitle(hb,'log \%','fontsize',labFont)
%titlemine('Weighted Score');

fname='cjf_rand_pref_weight_2dhist';
if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end


%% 
% %% plot 2D histogram of difference 
% hf=myFigure;
% imagesc(xl2,yl2,log10(bbird2./sum(sum(bbird2)).*100))
% hold on
% %plot([-0.5 21]./20,[-0.5 21]./20,'--k','linewidth',1.5)
% %plot([15.5 15.5]./20,[-0.5 21]./20,':k','linewidth',1.5)
% %plot([-0.5 21]./20,[15.5 15.5]./20,':k','linewidth',1.5)
% xlabelmine('Random',labFont);
% ylabelmine('Optimized-Random',labFont);
% %colormap(brewermap(256,'Reds'))
% colormap(flipud(magma(256)));
% %caxis(log10([0.01 20])); 
% % text(0.05,0.95,[num2str(q1,2) '\%'],'Interpreter','latex',...
% %     'fontsize',22)
% % text(0.8,0.95,[num2str(q2,2) '\%'],'Interpreter','latex',...
% %     'fontsize',22)
% % text(0.05,0.7,[num2str(q3,2) '\%'],'Interpreter','latex',...
% %     'fontsize',22)
% % text(0.8,0.7,[num2str(q4,2) '\%'],'Interpreter','latex',...
% %     'fontsize',22)
% 
% set(gca,'ydir','normal','fontsize',axFont)
% hb=colorbar;
% barTitle(hb,'log \%','fontsize',labFont)
% %titlemine('Weighted Score');
% 
% fname='cjf_rand_pref_weight_2dhist_diffR';
% if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end
% 
% %% plot 2D histogram of difference 
% hf=myFigure;
% imagesc(xl3,yl3,log10(bbird3./sum(sum(bbird3)).*100))
% hold on
% %plot([-0.5 21]./20,[-0.5 21]./20,'--k','linewidth',1.5)
% %plot([15.5 15.5]./20,[-0.5 21]./20,':k','linewidth',1.5)
% %plot([-0.5 21]./20,[15.5 15.5]./20,':k','linewidth',1.5)
% xlabelmine('Optimized',labFont);
% ylabelmine('Optimized-Random',labFont);
% %colormap(brewermap(256,'Reds'))
% colormap(flipud(magma(256)));
% %caxis(log10([0.01 20])); 
% % text(0.05,0.95,[num2str(q1,2) '\%'],'Interpreter','latex',...
% %     'fontsize',22)
% % text(0.8,0.95,[num2str(q2,2) '\%'],'Interpreter','latex',...
% %     'fontsize',22)
% % text(0.05,0.7,[num2str(q3,2) '\%'],'Interpreter','latex',...
% %     'fontsize',22)
% % text(0.8,0.7,[num2str(q4,2) '\%'],'Interpreter','latex',...
% %     'fontsize',22)
% 
% set(gca,'ydir','normal','fontsize',axFont)
% hb=colorbar;
% barTitle(hb,'log \%','fontsize',labFont)
% %titlemine('Weighted Score');
% 
% fname='cjf_rand_pref_weight_2dhist_diffP';
% if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end

%% Some stats 


ng=height(objectTableComp);
fprintf('Number of objects in random/preferred comparison: %i \n',ng)
fprintf('No. of objects with higher scores in preferred : %i , %3.1f %% \n',...
    sum(scP>scR),sum(scP>scR)/ng*100)
fprintf('No. of objects with lower scores in preferred: %i , %3.1f %%\n',...
    sum(scP<scR),sum(scP<scR)/ng*100)
fprintf('No. of objects with equal scores in preferred: %i , %3.1f %%\n',...
    sum(scP==scR),sum(scP==scR)/ng*100)
fprintf('No. of objects which are JF in preferred: %i , %3.1f %%\n',...
    sum(scP>=0.8),sum(scP>=0.8 )/ng*100)
fprintf('No. of objects which are JF in random : %i , %3.1f %%\n',...
    sum(scR>=0.8),sum( scR>=0.8)/ng*100)
fprintf('increase is  : %3.1f %%\n',...
    (sum(scP>=0.8)./sum(scR>=0.8)-1).*100)
fprintf('No. of objects which are JF in preferred but not in random: %i , %3.1f %%\n',...
    sum(scP>=0.8 & scR<0.8),sum(scP>=0.8 & scR<0.8)/ng*100)
fprintf('No. of objects which are JF in random but not in preferred: %i , %3.1f %%\n',...
    sum(scP<0.8 & scR>=0.8),sum(scP<0.8 & scR>=0.8)/ng*100)

