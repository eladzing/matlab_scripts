%% demographics of jellyfish by hosts
%% load data

global DEFAULT_MATFILE_DIR
load([ DEFAULT_MATFILE_DIR '\jf_statsByHosts_CJF.mat']);

global DEFAULT_PRINTOUT_DIR
outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/paper'];
%%
cmap=brewermap(256,'YlOrRd');
cmap(1,:)=[1 1 1];
cmap2=(viridis(256));
cmap2(1,:)=[1 1 1];
%%  
% figure
% semilogx(jfStats.M200c,jfStats.sampleSats,'.');
% xlabelmine('host mass');
% ylabelmine('No. of Sats');
% 
% figure
% loglog(jfStats.M200c,jfStats.sampleSats,'.');
% xlabelmine('host mass');
% ylabelmine('No. of Sats');

%% build birds 
[bird, binsize, xl,yl]= histogram2d(log10(jfStats.M200c),jfStats.sampleSats,ones(size(jfStats.sampleSats)),...
    'xlim',[10 15],'ylim',[0 65],'len',65);
mm=jfStats.sim=="TNG100";
[bird100, ~, ~,~]= histogram2d(log10(jfStats.M200c(mm)),jfStats.sampleSats(mm),ones(size(jfStats.sampleSats(mm))),...
    'xlim',[10 15],'ylim',[0 65],'len',65);
[bird50, ~, ~,~]= histogram2d(log10(jfStats.M200c(~mm)),jfStats.sampleSats(~mm),ones(size(jfStats.sampleSats(~mm))),...
    'xlim',[10 15],'ylim',[0 65],'len',65);
%%

% myFigure;
% %t=tiledlayout(6,5);
% 
% %nexttile([6 3])
% imagesc(xl,yl,log10(squeeze(bird(:,:,1))))
% set(gca,'ydir','normal','fontsize',14)
% xlabelmine('log Host Mass');
% ylabelmine('No. of Satellites');
% colormap(cmap2)
% caxis(log10([0.9 1000]))
% cbh=colorbar;
cax=log10( [0.8 1000]);

myFigure;
%nexttile([3 2])
countMap=squeeze(bird100(:,:,1));
imagesc(xl,yl,log10(countMap));
set(gca,'ydir','normal','fontsize',14)
text(10.2,60,'TNG100','Fontsize',16,'interpreter','latex')
text(10.2,55,num2str(sum(sum(countMap))),'Fontsize',16,'interpreter','latex')
xlabelmine('log Host Mass');
ylabelmine('No. of Sample Satellites');
colormap(cmap2)
caxis(cax);
cbh(1)=colorbar('fontsize',14);
barTitle(cbh(1),'log No. of Hosts')
fname='cjf_satInHost_TNG100';
printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir);
%titlemine('TNG100',12);


myFigure;
%nexttile([3 2])
countMap=squeeze(bird50(:,:,1));
imagesc(xl,yl,log10(countMap));
set(gca,'ydir','normal','fontsize',14)
text(10.2,60,'TNG50','Fontsize',16,'interpreter','latex')
text(10.2,55,num2str(sum(sum(countMap))),'Fontsize',16,'interpreter','latex')
xlabelmine('log Host Mass');
ylabelmine('No. of Sample Satellites');
colormap(cmap2)
caxis(cax);
cbh=colorbar('fontsize',14);
barTitle(cbh,'log No. of Hosts')
fname='cjf_satInHost_TNG50';
printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir);
%titlemine('TNG100',12);

%titlemine('TNG50',12);
%cbh.Layout.Tile = 'east';
%%
zr=illustris.utils.snap2redshift(jfStats.snap);
gr=ones(size(zr));
gr(zr<=0.5)=0;
gr(zr>0.5 & zr<=1)=1;
gr(zr>1 & zr<=2.1)=2;

for j=0:2
mm=jfStats.sim=="TNG100" & gr==j;
[bird, ~, ~,~]= histogram2d(log10(jfStats.M200c(mm)),jfStats.sampleSats(mm),ones(size(jfStats.sampleSats(mm))),...
    'xlim',[10 15],'ylim',[0 65],'len',65);
% [bird50, ~, ~,~]= histogram2d(log10(jfStats.M200c(~mm)),jfStats.sampleSats(~mm),ones(size(jfStats.sampleSats(~mm))),...
%     'xlim',[10 15],'ylim',[0 65],'len',65);

cntMap(:,:,j+1)=squeeze(bird(:,:,1));
end

%%
myFigure('pos',[991   179   875   720]);
t=tiledlayout(6,3);
nexttile([4 3])
countMap=squeeze(bird100(:,:,1));
imagesc(xl,yl,log10(countMap));
set(gca,'ydir','normal','fontsize',16)
text(10.2,60,'TNG100','Fontsize',18,'interpreter','latex')
text(10.2,55,num2str(sum(sum(countMap))),'Fontsize',18,'interpreter','latex')
%xlabelmine('log Host Mass');
%ylabelmine('No. of Sample Satellites');
colormap(cmap2)
caxis(cax);
cbh=colorbar('fontsize',16);
barTitle(cbh(1),'$\log N$','fontsize',18)
%tag={'$z=[0\; 0.5]$','$z=[0.5\; 1]$','$z=[1 2]$'};
tag={'$0\leq z \leq 0.5$','$0.5< z \leq 1$','$1< z \leq 2$'};
for k=1:3
    
nexttile([2,1])
ccm=squeeze(cntMap(:,:,k));
imagesc(xl,yl,log10(ccm));
set(gca,'ydir','normal','fontsize',12)
text(10.2,55,tag{k},'Fontsize',16,'interpreter','latex')
text(10.2,40,num2str(sum(sum(ccm))),'Fontsize',16,'interpreter','latex')
colormap(cmap2)
caxis(cax);
end
t.XLabel.String='log Host Mass';
t.XLabel.FontSize=18;
t.XLabel.Interpreter='latex';
t.YLabel.String='No. of Sample Satellites';
t.YLabel.FontSize=18;
t.YLabel.Interpreter='latex';
t.TileSpacing='compact';


cbh.Layout.Tile = 'east';
%t.YLabel('No. of Sample Satellites');


fname='cjf_satInHost_2_TNG100';
printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir);
%titlemine('TNG100',12);


%% TNG50 

for j=0:2
mm=jfStats.sim=="TNG50" & gr==j;
[bird, ~, ~,~]= histogram2d(log10(jfStats.M200c(mm)),jfStats.sampleSats(mm),ones(size(jfStats.sampleSats(mm))),...
    'xlim',[10 15],'ylim',[0 65],'len',65);
% [bird50, ~, ~,~]= histogram2d(log10(jfStats.M200c(~mm)),jfStats.sampleSats(~mm),ones(size(jfStats.sampleSats(~mm))),...
%     'xlim',[10 15],'ylim',[0 65],'len',65);

cntMap(:,:,j+1)=squeeze(bird(:,:,1));
end

%%
myFigure('pos',[991   179   875   720]);
t=tiledlayout(6,3);
nexttile([4 3])
countMap=squeeze(bird50(:,:,1));
imagesc(xl,yl,log10(countMap));
set(gca,'ydir','normal','fontsize',16)
text(10.2,60,'TNG50','Fontsize',18,'interpreter','latex')
text(10.2,55,num2str(sum(sum(countMap))),'Fontsize',18,'interpreter','latex')
%xlabelmine('log Host Mass');
%ylabelmine('No. of Sample Satellites');
colormap(cmap2)
caxis(cax);
cbh=colorbar('fontsize',16);
barTitle(cbh(1),'$\log N$','fontsize',18)
%tag={'$z=[0\; 0.5]$','$z=[0.5\; 1]$','$z=[1 2]$'};
tag={'$0\leq z \leq 0.5$','$0.5< z \leq 1$','$1< z \leq 2$'};
for k=1:3
    
nexttile([2,1])
ccm=squeeze(cntMap(:,:,k));
imagesc(xl,yl,log10(ccm));
set(gca,'ydir','normal','fontsize',11)
text(10.2,55,tag{k},'Fontsize',16,'interpreter','latex')
text(10.2,40,num2str(sum(sum(ccm))),'Fontsize',16,'interpreter','latex')
colormap(cmap2)
caxis(cax);
end
t.XLabel.String='log Host Mass';
t.XLabel.FontSize=18;
t.XLabel.Interpreter='latex';
t.YLabel.String='No. of Sample Satellites';
t.YLabel.FontSize=18;
t.YLabel.Interpreter='latex';
t.TileSpacing='compact';


cbh.Layout.Tile = 'east';
%t.YLabel('No. of Sample Satellites');


fname='cjf_satInHost_2_TNG50';
printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir);

 %%

% figure
% scatterhist(log10(jfStats.M200c),log10(jfStats.sampleSats),...
%     'group',jfStats.sim=="TNG100",...
%     'location','northeast','Direction','out','kernel','off',...
%     'nbins',30);
% xlabelmine('host mass');
% ylabelmine('No. of Sats');

%%
% zr=illustris.utils.snap2redshift(jfStats.snap);
% gr=ones(size(zr));
% gr(zr<=0.5)=0;
% gr(zr>0.5 & zr<=1)=1;
% gr(zr>1 & zr<=2)=2;
% 
% for j=1:3
%     switch j
%         case 1
%             msk=true(size(jfStats.sim));
%             tit='All';
%         case 2
%             msk=jfStats.sim=="TNG100";
%             tit='TNG100';
%         case 3
%             msk=jfStats.sim=="TNG50";
%             tit='TNG50';
%             
%     end
%     
%     figure
%     scatterhist(log10(jfStats.M200c(msk)),log10(jfStats.sampleSats(msk)),...
%         'group',gr(msk),...
%         'location','northeast','Direction','out','kernel','off',...
%         'nbins',30,'markersize',5);
%     xlabelmine('host mass');
%     ylabelmine('No. of Sats');
%     titlemine(tit);
% end

%% plot jf fractions

zr=illustris.utils.snap2redshift(jfStats.snap);
gr=ones(size(zr));
gr(zr<=0.5)=0;
gr(zr>0.5 & zr<=1)=1;
gr(zr>1 & zr<=2.1)=2;
% grr=string(gr);
% grr(zr<=0.5)="$z\leq 0.5$";
% grr(zr>0.5 & zr<=1)="$z>0.5 \& z\leq 1$";
% grr(zr>1 & zr<=2.1)="$z>1 \& z\leq 2$";

jff=jfStats.JFNum./jfStats.sampleSats;
mmm=jff==0;

% msk=~mmm;
% tab=table(jfStats.M200c(msk),jff(msk),grr(msk),'VariableNames',{"host mass","JF Fraction","redshift"});
%
%  msk=jfStats.sim=="TNG100" & ~mmm;
% tab100=table(jfStats.M200c(msk),jff(msk),grr(msk),'VariableNames',{'host mass','JF Fraction','redshift'});
%
% msk=jfStats.sim=="TNG50" & ~mmm ;
% tab50=table(jfStats.M200c(msk),jff(msk),grr(msk),'VariableNames',{'host mass','JF Fraction','redshift'});



% give fictitios frac for zero
% base=1e-2;
% scat=0.2;
% mmm=jff==0;
% jff(mmm)=base.*10.^(scat.*rand(1,sum(mmm)));

cmap=brewermap(5,'Dark2');


msk100=jfStats.sim=="TNG100"  & ~mmm;
msk50=jfStats.sim=="TNG50" & ~mmm ;
%%

hf=myFigure ;
  
   scatterhist(log10(jfStats.M200c(msk50)),log10(jff(msk50)),...
        'group',gr(msk50),'legend','off','style','stairs',...
        'location','Northwest','Direction','out','kernel','off',...
        'nbins',10,'markersize',[4,6,8],'marker','osd','color',cmap([5 2 3],:));
    set(gca,'xaxislocation','bottom',...
        'yaxislocation','right','fontsize',16)
    
    % set linwidth for histogram
    for i=2:3
        for j=1:3
            hf.Children(i).Children(j).LineWidth=1.8;
        end
    end
   
    
    %hf.Children(2).XScale='log';
    %hf.Children(3).XScale='log';
    
    hf.Children(5).Children(1).MarkerFaceColor=cmap(3,:);
    hf.Children(5).Children(2).MarkerFaceColor=cmap(2,:);
    hf.Children(5).Children(3).MarkerFaceColor=cmap(5,:);
    
    hf.Children(4).String = {'z=[0 0.5]','z=[0.5 1]','z=[1 2]'};  %{"$z\leq 0.5$","$z>0.5 \& z\leq 1$","$z>1 \& z\leq 2$"};
    hf.Children(4).Interpreter ='latex';
    hf.Children(4).FontSize =16 ;

  grid
    % hold on
    % plot(log10(jfStats.M200c(msk & mmm)),log10(jff(msk & mmm)),'.')
    text(11.6,-1.4,"TNG50",'Interpreter','latex','fontsize',16);
    text(11.6,-1.6,num2str(sum(msk50 & gr==0)),'Interpreter','latex','Color',cmap(5,:),'fontsize',16);
    text(11.6,-1.75,num2str(sum(msk50 & gr==1)),'Interpreter','latex','Color',cmap(2,:),'fontsize',16);
    text(11.6,-1.9,num2str(sum(msk50 & gr==2)),'Interpreter','latex','Color',cmap(3,:),'fontsize',16);
    xlim([11.4 14.7])
    ylim([-2 0.05])
    xlabelmine('log Host Mass',18);
    ylabelmine('JF Fraction',18);
    
    fname='cjf_jfhost_jfFrac_TNG50';
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir);
    
    
    %%
    hf=myFigure ;
  
   scatterhist(log10(jfStats.M200c(msk100)),log10(jff(msk100)),...
        'group',gr(msk100),'legend','on','style','stairs',...
        'location','Northeast','Direction','out','kernel','off',...
        'nbins',10,'markersize',[4,6,8],'marker','osd','color',cmap([5 2 3],:));
    set(gca,'xaxislocation','bottom',...
        'yaxislocation','left','fontsize',16)
    
    % set linwidth for histogram
    for i=2:3
        for j=1:3
            hf.Children(i).Children(j).LineWidth=1.8;
        end
    end
   
    
    %hf.Children(2).XScale='log';
    %hf.Children(3).XScale='log';
    
    hf.Children(5).Children(1).MarkerFaceColor=cmap(3,:);
    hf.Children(5).Children(2).MarkerFaceColor=cmap(2,:);
    hf.Children(5).Children(3).MarkerFaceColor=cmap(5,:);
    
    hf.Children(4).String = {'$0\leq z \leq 0.5$','$0.5< z \leq 1$','$1< z \leq 2$'};  %{"$z\leq 0.5$","$z>0.5 \& z\leq 1$","$z>1 \& z\leq 2$"};
    hf.Children(4).Interpreter ='latex';
    hf.Children(4).FontSize =16 ;

  grid
    % hold on
    % plot(log10(jfStats.M200c(msk & mmm)),log10(jff(msk & mmm)),'.')
    text(11.6,-1.4,"TNG100",'Interpreter','latex','fontsize',16);
    text(11.6,-1.6,num2str(sum(msk100 & gr==0)),'Interpreter','latex','Color',cmap(5,:),'fontsize',16);
    text(11.6,-1.75,num2str(sum(msk100 & gr==1)),'Interpreter','latex','Color',cmap(2,:),'fontsize',16);
    text(11.6,-1.9,num2str(sum(msk100 & gr==2)),'Interpreter','latex','Color',cmap(3,:),'fontsize',16);
    xlim([11.4 14.7])
    ylim([-2 0.05])
    xlabelmine('log Host Mass',18);
    ylabelmine('',18);
    %ylabelmine('JF Fraction',18);
    
    fname='cjf_jfhost_jfFrac_TNG100';
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir);
    %%
    
    
    
    
   scatterhist(log10(jfStats.M200c(msk100)),log10(jff(msk100)),...
        'group',gr(msk100),'legend','off',...
        'location','Northeast','Direction','out','kernel','on',...
        'nbins',50,'markersize',4,'marker','.od','color',cmap([1 2 5],:));
    set(gca,'xaxislocation','bottom',...
        'yaxislocation','left')
   
   
%     scatterhist((jfStats.M200c(mskk)),(jff(mskk)),...
%         'group',gr(mskk),...
%         'location','northeast','Direction','out','kernel','off',...
%         'nbins',50,'markersize',4,'marker','x+d');
%     set(gca,'yscale','log','xscale','log','xaxislocation','bottom',...
%         'yaxislocation','left')
%  scatterhist(log10(jfStats.M200c(mskk)),log10(jff(mskk)),...
%         'group',gr(mskk),'legend',lflag,...
%         'location',loc,'Direction','out','kernel','on',...
%         'nbins',50,'markersize',4,'marker','xod','color',cmap([1 2 5],:));
%     set(gca,'xaxislocation','bottom',...
%         'yaxislocation',aloc)
    colormap(cmap)
    %hf.Children(2).XScale='log';
    %hf.Children(3).XScale='log';
    %hf.Children(4).String = {'z=[0 0.5]','z=[0.5 1]','z=[1 2]'};  %{"$z\leq 0.5$","$z>0.5 \& z\leq 1$","$z>1 \& z\leq 2$"};
    
    
    %
    %hf.Children(4).Interpreter ='latex';
    %hf.Children(4).FontSize =16 ;
    
    grid
    % hold on
    % plot(log10(jfStats.M200c(msk & mmm)),log10(jff(msk & mmm)),'.')
    text(10.6,-0.2,tit,'Interpreter','latex');
    text(10.6,-0.4,num2str(sum(mm0)),'Interpreter','latex','Color',cmap(1,:));
    text(10.6,-0.55,num2str(sum(mm1)),'Interpreter','latex','Color',cmap(2,:));
    text(10.6,-0.7,num2str(sum(mm2)),'Interpreter','latex','Color',cmap(5,:));
    xlim([10.3 14.7])
    ylim([-2 0.05])
    xlabelmine('host mass');
    ylabelmine('JF Fraction');
    
    
    
    
    
    

%%

figure
histogram(log10(jfStats.M200c(mmm)),linspace(10.4,15,50),'facecolor',cmap(2,:));
hold on;
histogram(log10(jfStats.M200c(~mmm)),linspace(10.4,15,50),'facecolor',cmap(1,:));
%text(10^14.2,1000,num2str(sum(mmm)),'Interpreter','latex','Color',cmap(2,:));
%text(10^14.2,700,num2str(sum(~mmm)),'Interpreter','latex','Color',cmap(1,:));
legend({['No JF (' num2str(sum(mmm)) ')'], ['JF (' num2str(sum(~mmm)) ')'] },'Interpreter','latex','fontsize',16    )
set(gca,'yscale','log')
ylim([0.5 2000 ])

%%

mmm=jff==0;
m100=jfStats.sim=='TNG100';
m50=jfStats.sim=='TNG50';


myFigure;
histogram(log10(jfStats.M200c(mmm & m100)),linspace(10.4,15,50),'facecolor',cmap(2,:));
hold on;
histogram(log10(jfStats.M200c(~mmm & m100)),linspace(10.4,15,50),'facecolor',cmap(1,:));
histogram(log10(jfStats.M200c(mmm)),linspace(10.4,15,50),'facecolor','none',...
    'edgecolor',cmap(2,:),'linewidth',1.1);
histogram(log10(jfStats.M200c(~mmm)),linspace(10.4,15,50),'facecolor','none',...
    'edgecolor',cmap(1,:),'linewidth',1.1);
%text(10^14.2,1000,num2str(sum(mmm)),'Interpreter','latex','Color',cmap(2,:));
%text(10^14.2,700,num2str(sum(~mmm)),'Interpreter','latex','Color',cmap(1,:));
legend({['No JF (' num2str(sum(mmm & m100)) ')'], ['JF (' num2str(sum(~mmm & m100)) ')'] },'Interpreter','latex','fontsize',16    )
xlabelmine('log Host Mass');
ylabelmine('No. of Hosts');
set(gca,'yscale','log','fontsize',14)
titlemine('TNG100',16)
%text(14.3,400,'TNG100','Fontsize',16,'interpreter','latex')
ylim([0.5 2000 ])
fname='cjf_jfhost_hist_TNG100';
printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir);

myFigure;
histogram(log10(jfStats.M200c(mmm & m50)),linspace(10.4,15,50),'facecolor',cmap(2,:));
hold on;
histogram(log10(jfStats.M200c(~mmm & m50)),linspace(10.4,15,50),'facecolor',cmap(1,:));

histogram(log10(jfStats.M200c(mmm)),linspace(10.4,15,50),'facecolor','none',...
    'edgecolor',cmap(2,:),'linewidth',1.1);
histogram(log10(jfStats.M200c(~mmm)),linspace(10.4,15,50),'facecolor','none',...
    'edgecolor',cmap(1,:),'linewidth',1.1);
%text(10^14.2,1000,num2str(sum(mmm)),'Interpreter','latex','Color',cmap(2,:));
%text(10^14.2,700,num2str(sum(~mmm)),'Interpreter','latex','Color',cmap(1,:));
legend({['No JF (' num2str(sum(mmm & m50)) ')'], ['JF (' num2str(sum(~mmm & m50)) ')'] },'Interpreter','latex','fontsize',16    )
xlabelmine('log Host Mass');
ylabelmine('No. of Hosts');
set(gca,'yscale','log','fontsize',14)
%text(14.3,400,'TNG50','Fontsize',16,'interpreter','latex')
titlemine('TNG50',16)
ylim([0.5 2000 ])
fname='cjf_jfhost_hist_TNG50';
printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir);