%% demographics of jellyfish by hosts
%% load data
if setupFlag
    global DEFAULT_MATFILE_DIR
    load([ DEFAULT_MATFILE_DIR '\jf_statsByHosts_CJF.mat']);
    
    global DEFAULT_PRINTOUT_DIR
    outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/paper'];
end
%%
cmap=brewermap(256,'YlOrRd');
cmap(1,:)=[1 1 1];
cmap2=(viridis(256));
cmap2(1,:)=[1 1 1];

 axFont=18;
 legFont=20;
 labFont=20;
 
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
% xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$);
% ylabelmine('No. of Satellites');
% colormap(cmap2)
% caxis(log10([0.9 1000]))
% cbh=colorbar;
cax=log10( [0.8 1000]);

myFigure;
%nexttile([3 2])
countMap=squeeze(bird100(:,:,1));
imagesc(xl,yl,log10(countMap));
set(gca,'ydir','normal','fontsize',axFont)
text(10.2,60,'TNG100','Fontsize',16,'interpreter','latex')
text(10.2,55,num2str(sum(sum(countMap))),'Fontsize',16,'interpreter','latex')
xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$');
ylabelmine('No. of Inspected Satellites');
colormap(cmap2)
caxis(cax);
cbh(1)=colorbar('fontsize',axFont);
barTitle(cbh(1),'log No. of Hosts')
fname='cjf_satInHost_TNG100';
if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end
%titlemine('TNG100',12);


myFigure;
%nexttile([3 2])
countMap=squeeze(bird50(:,:,1));
imagesc(xl,yl,log10(countMap));
set(gca,'ydir','normal','fontsize',axFont)
text(10.2,60,'TNG50','Fontsize',16,'interpreter','latex')
text(10.2,55,num2str(sum(sum(countMap))),'Fontsize',16,'interpreter','latex')
xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$');
ylabelmine('No. of Inspected Satellites');
colormap(cmap2)
caxis(cax);
cbh=colorbar('fontsize',axFont);
barTitle(cbh,'log No. of Hosts')
fname='cjf_satInHost_TNG50';
if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end
%titlemine('TNG100',12);

%titlemine('TNG50',12);
%cbh.Layout.Tile = 'east';
%%
zr=round(illustris.utils.snap2redshift(jfStats.snap),2);
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
set(gca,'ydir','normal','fontsize',axFont)
text(10.2,60,'TNG100','Fontsize',legFont,'interpreter','latex')
text(10.2,55,num2str(sum(sum(countMap))),'Fontsize',legFont,'interpreter','latex')
%xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$');
%ylabelmine('No. of Inspected Satellites');
colormap(cmap2)
caxis(cax);
cbh=colorbar('fontsize',axFont);
barTitle(cbh(1),'$\log N$','fontsize',labFont)
%tag={'$z=[0\; 0.5]$','$z=[0.5\; 1]$','$z=[1 2]$'};
tag={'$0\leq z \leq 0.5$','$0.5< z \leq 1$','$1< z \leq 2$'};
for k=1:3
    
nexttile([2,1])
ccm=squeeze(cntMap(:,:,k));
imagesc(xl,yl,log10(ccm));
set(gca,'ydir','normal','fontsize',12)
text(10.2,55,tag{k},'Fontsize',18,'interpreter','latex')
text(10.2,40,num2str(sum(sum(ccm))),'Fontsize',18,'interpreter','latex')
colormap(cmap2)
caxis(cax);
end
t.XLabel.String='log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$';
t.XLabel.FontSize=labFont;
t.XLabel.Interpreter='latex';
t.YLabel.String='No. of Inspected Satellites';
t.YLabel.FontSize=labFont;
t.YLabel.Interpreter='latex';
t.TileSpacing='compact';


cbh.Layout.Tile = 'east';
%t.YLabel('No. of Inspected Satellites');


fname='cjf_satInHost_2_TNG100';
if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end
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
set(gca,'ydir','normal','fontsize',axFont)
text(10.2,60,'TNG50','Fontsize',legFont,'interpreter','latex')
text(10.2,55,num2str(sum(sum(countMap))),'Fontsize',legFont,'interpreter','latex')
%xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$');
%ylabelmine('No. of Inspected Satellites');
colormap(cmap2)
caxis(cax);
cbh=colorbar('fontsize',axFont);
barTitle(cbh(1),'$\log N$','fontsize',axFont)
%tag={'$z=[0\; 0.5]$','$z=[0.5\; 1]$','$z=[1 2]$'};
tag={'$0\leq z \leq 0.5$','$0.5< z \leq 1$','$1< z \leq 2$'};
for k=1:3
    
nexttile([2,1])
ccm=squeeze(cntMap(:,:,k));
imagesc(xl,yl,log10(ccm));
set(gca,'ydir','normal','fontsize',12)
text(10.2,55,tag{k},'Fontsize',18,'interpreter','latex')
text(10.2,40,num2str(sum(sum(ccm))),'Fontsize',18,'interpreter','latex')
colormap(cmap2)
caxis(cax);
end
t.XLabel.String='log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$';
t.XLabel.FontSize=labFont;
t.XLabel.Interpreter='latex';
t.YLabel.String='No. of Inspected Satellites';
t.YLabel.FontSize=labFont;
t.YLabel.Interpreter='latex';
t.TileSpacing='compact';


cbh.Layout.Tile = 'east';
%t.YLabel('No. of Inspected Satellites');


fname='cjf_satInHost_2_TNG50';
if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end

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

zr=round(illustris.utils.snap2redshift(jfStats.snap),2);
gr=ones(size(zr));
gr(zr<=0.5)=0;
gr(zr>0.5 & zr<=1)=1;
gr(zr>1 & zr<=2.1)=2;
% grr=string(gr);
% grr(zr<=0.5)="$z\leq 0.5$";
% grr(zr>0.5 & zr<=1)="$z>0.5 \& z\leq 1$";
% grr(zr>1 & zr<=2.1)="$z>1 \& z\leq 2$";

jff=jfStats.JFNumWeighted./jfStats.sampleSats;
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

cmap=brewermap(12,'Paired');
%cmap=brewermap(5,'Dark2');


msk100=jfStats.sim=="TNG100"  & ~mmm;
msk50=jfStats.sim=="TNG50" & ~mmm ;

%% print out some stats on JF fractions 
qn50=[];
qn100=[];
for i=[0 1 2]
    mask50=gr==i & msk50;
    mask100=gr==i & msk100;
    qn50(i+1,:)=quantile(jff(mask50),[0.05 0.1 0.2 0.25 0.5 0.75 0.8 0.9 0.95]);
    qn100(i+1,:)=quantile(jff(mask100),[0.05 0.1 0.2 0.25 0.5 0.75 0.8 0.9 0.95]);
end

%myFigure
figure(2)
xx=[0  0.5 1 2];

 stairs(xx,qn50([1  2 3 3],4),'-r','linewidth',2.5)
 hold on
 stairs(xx,qn50([1  2 3 3],2),'-r','linewidth',0.9)
 stairs(xx,qn50([1  2 3 3],6),'-r','linewidth',0.9)
% 
% stairs(xx,qn100([1  2 3 3],4),'-b','linewidth',2.5)
% stairs(xx,qn100([1  2 3 3],2),'-b','linewidth',0.9)
% stairs(xx,qn100([1  2 3 3],6),'-b','linewidth',0.9)
% xlim([0 2]);
% ylim([0 0.35])
% xlabelmine('redshift');
% ylabelmine('median JF fraction');
%%

grr=ones(size(zr))*-1;
grr(zr<=0.5 & msk50)=0;
grr(zr<=0.5 & msk100)=1;
grr((zr>0.5 & zr<=1) & msk50)=10 ;
grr((zr>0.5 & zr<=1) & msk100)=11 ;
grr((zr>1 & zr<=2.1) & msk50)=20;
grr((zr>1 & zr<=2.1) & msk100)=21;
% grr=string(gr);
% grr(zr<=0.5)="$z\leq 0.5$";
% grr(zr>0.5 & zr<=1)="$z>0.5 \& z\leq 1$";
% grr(zr>1 & zr<=2.1)="$z>1 \& z\leq 2$";
colset1=brewermap(2,'Set1');
cols(1,:)=colset1(1,:);
cols(3,:)=colset1(1,:);
cols(5,:)=colset1(1,:);
cols(2,:)=colset1(2,:);
cols(4,:)=colset1(2,:);
cols(6,:)=colset1(2,:);

%%

myFigure;%('pos',[991   179   875   720]);
hb=boxplot(jff(grr~=-1),grr(grr~=-1),'PlotStyle','traditional','colors','rbrbrb','jitter',1.2,...
    'symbol','x','colors',cols,...
    'labels',{'',' ',' ','','',' '});
set(hb,'linewidth',1.5,'markersize',8)
grid
text(5.2,0.9,'TNG50','color',colset1(1,:),'fontsize',legFont,'Interpreter','latex')
text(5.2,0.81,'TNG100','color',colset1(2,:),'fontsize',legFont,'Interpreter','latex')
text(1.1,-0.1,'$z\le0.5$','Interpreter','latex','fontsize',18);
text(2.9,-0.1,'$0.5<z\le 1$','Interpreter','latex','fontsize',18);
text(5,-0.1,'$1< z\le 2$','Interpreter','latex','fontsize',18);
%boxplot(jff(grr~=-1),grr(grr~=-1),'whisker',1,'PlotStyle','compact');
set(gca,'fontsize',16)
%xlabelmine('Redshift Bins')
ylabelmine('JF Fraction',labFont);


  fname='cjf_jfhost_jfFrac_boxplot_zredBin';
    if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end
    

   % fprintf('jf fraction quantilies over hosts with jf: %s \n',num2str(round(median(jff(mask50)),2)))
% 
%%

hf=myFigure('pos',[97   444   838   687]);
  
   scatterhist(log10(jfStats.M200c(msk50)),log10(jff(msk50)),...
        'group',gr(msk50),'legend','off','style','stairs',...
        'location','Northwest','Direction','out','kernel','off',...
        'nbins',10,'markersize',[5,8,10],'marker','osd','color',cmap([4 8 10 ],:));%[5 2 3]
    set(gca,'xaxislocation','bottom',...
        'yaxislocation','right','fontsize',axFont)
    
    % set linwidth for histogram
    for i=2:3
        for j=1:3
            hf.Children(i).Children(j).LineWidth=1.8;
        end
    end
   
    
    %hf.Children(2).XScale='log';
    %hf.Children(3).XScale='log';
    
    hf.Children(5).Children(1).MarkerFaceColor=cmap(10  ,:);
    hf.Children(5).Children(2).MarkerFaceColor=cmap(8,:);
    hf.Children(5).Children(3).MarkerFaceColor=cmap(4,:);
    
    hf.Children(4).String = {'z=[0 0.5]','z=[0.5 1]','z=[1 2]'};  %{"$z\leq 0.5$","$z>0.5 \& z\leq 1$","$z>1 \& z\leq 2$"};
    hf.Children(4).Interpreter ='latex';
    hf.Children(4).FontSize =legFont ;

  grid
    % hold on
    % plot(log10(jfStats.M200c(msk & mmm)),log10(jff(msk & mmm)),'.')
    text(10.6,-1.4,"TNG50",'Interpreter','latex','fontsize',legFont);
    text(10.6,-1.6,num2str(sum(msk50 & gr==0)),'Interpreter','latex','Color',cmap(4,:),'fontsize',legFont);
    text(10.6,-1.75,num2str(sum(msk50 & gr==1)),'Interpreter','latex','Color',cmap(8,:),'fontsize',legFont);
    text(10.6,-1.9,num2str(sum(msk50 & gr==2)),'Interpreter','latex','Color',cmap(10,:),'fontsize',legFont);
    xlim([10.4 14.5])
    %xlim([11.4 14.7])
    ylim([-2 0.05])
    xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$',labFont);
    ylabelmine('log JF Fraction',labFont);
    
    fname='cjf_jfhost_jfFrac_TNG50';
    if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end
    
    
    %%
    hf=myFigure('pos',[97   444   838   687]);
  
   scatterhist(log10(jfStats.M200c(msk100)),log10(jff(msk100)),...
        'group',gr(msk100),'legend','on','style','stairs',...
        'location','Northeast','Direction','out','kernel','off',...
        'nbins',10,'markersize',[5,8,10],'marker','osd','color',cmap([4 8 10],:));
    set(gca,'xaxislocation','bottom',...
        'yaxislocation','left','fontsize',axFont)
    
    % set linwidth for histogram
    for i=2:3
        for j=1:3
            hf.Children(i).Children(j).LineWidth=1.8;
        end
    end
   
    
    %hf.Children(2).XScale='log';
    %hf.Children(3).XScale='log';
    
    hf.Children(5).Children(1).MarkerFaceColor=cmap(10,:);
    hf.Children(5).Children(2).MarkerFaceColor=cmap(8,:);
    hf.Children(5).Children(3).MarkerFaceColor=cmap(4,:);
    
    hf.Children(4).String = {'$0\leq z \leq 0.5$','$0.5< z \leq 1$','$1< z \leq 2$'};  %{"$z\leq 0.5$","$z>0.5 \& z\leq 1$","$z>1 \& z\leq 2$"};
    hf.Children(4).Interpreter ='latex';
    hf.Children(4).FontSize =legFont ;
    hf.Children(4).Box='off';
  grid
    % hold on
    % plot(log10(jfStats.M200c(msk & mmm)),log10(jff(msk & mmm)),'.')
    text(12.18,-1.4,"TNG100",'Interpreter','latex','fontsize',legFont);
    text(12.18,-1.6,num2str(sum(msk100 & gr==0)),'Interpreter','latex','Color',cmap(4,:),'fontsize',legFont);
    text(12.18,-1.75,num2str(sum(msk100 & gr==1)),'Interpreter','latex','Color',cmap(8,:),'fontsize',legFont);
    text(12.18,-1.9,num2str(sum(msk100 & gr==2)),'Interpreter','latex','Color',cmap(10,:),'fontsize',legFont);
    xlim([12 15])
    %xlim([11.4 14.7])
    ylim([-2 0.05])
    xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$',labFont);
    ylabelmine(' ',labFont);
    %ylabelmine('JF Fraction',18);
    
    fname='cjf_jfhost_jfFrac_TNG100';
    if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end
    %%
    
    
%     
%     
%    scatterhist(log10(jfStats.M200c(msk100)),log10(jff(msk100)),...
%         'group',gr(msk100),'legend','off',...
%         'location','Northeast','Direction','out','kernel','on',...
%         'nbins',50,'markersize',4,'marker','.od','color',cmap([1 2 5],:));
%     set(gca,'xaxislocation','bottom',...
%         'yaxislocation','left')
%    
%    
% %     scatterhist((jfStats.M200c(mskk)),(jff(mskk)),...
% %         'group',gr(mskk),...
% %         'location','northeast','Direction','out','kernel','off',...
% %         'nbins',50,'markersize',4,'marker','x+d');
% %     set(gca,'yscale','log','xscale','log','xaxislocation','bottom',...
% %         'yaxislocation','left')
% %  scatterhist(log10(jfStats.M200c(mskk)),log10(jff(mskk)),...
% %         'group',gr(mskk),'legend',lflag,...
% %         'location',loc,'Direction','out','kernel','on',...
% %         'nbins',50,'markersize',4,'marker','xod','color',cmap([1 2 5],:));
% %     set(gca,'xaxislocation','bottom',...
% %         'yaxislocation',aloc)
%     colormap(cmap)
%     %hf.Children(2).XScale='log';
%     %hf.Children(3).XScale='log';
%     %hf.Children(4).String = {'z=[0 0.5]','z=[0.5 1]','z=[1 2]'};  %{"$z\leq 0.5$","$z>0.5 \& z\leq 1$","$z>1 \& z\leq 2$"};
%     
%     
%     %
%     %hf.Children(4).Interpreter ='latex';
%     %hf.Children(4).FontSize =16 ;
%     
%     grid
%     % hold on
%     % plot(log10(jfStats.M200c(msk & mmm)),log10(jff(msk & mmm)),'.')
%     text(10.6,-0.2,tit,'Interpreter','latex');
%     text(10.6,-0.4,num2str(sum(mm0)),'Interpreter','latex','Color',cmap(1,:));
%     text(10.6,-0.55,num2str(sum(mm1)),'Interpreter','latex','Color',cmap(2,:));
%     text(10.6,-0.7,num2str(sum(mm2)),'Interpreter','latex','Color',cmap(5,:));
%     xlim([10.3 14.7])
%     ylim([-2 0.05])
%     xlabelmine('host mass');
%     ylabelmine('log JF Fraction');
%     
    
    
    
    
    

%%

myFigure;
histogram(log10(jfStats.M200c(mmm)),linspace(10.4,15,50),'facecolor',cmap(2,:));
hold on;
histogram(log10(jfStats.M200c(~mmm)),linspace(10.4,15,50),'facecolor',cmap(6,:));
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
histogram(log10(jfStats.M200c(mmm & m100)),linspace(10.4,15,50),'facecolor',cmap(6,:));
hold on;
histogram(log10(jfStats.M200c(~mmm & m100)),linspace(10.4,15,50),'facecolor',cmap(2,:));
histogram(log10(jfStats.M200c(mmm)),linspace(10.4,15,50),'displayStyle','stairs',...
    'edgecolor',cmap(6,:),'linewidth',1.5);
histogram(log10(jfStats.M200c(~mmm)),linspace(10.4,15,50),'displayStyle','stairs',...
    'edgecolor',cmap(2,:),'linewidth',1.5);
%text(10^14.2,1000,num2str(sum(mmm)),'Interpreter','latex','Color',cmap(2,:));
%text(10^14.2,700,num2str(sum(~mmm)),'Interpreter','latex','Color',cmap(1,:));
legend({['Non-JF (' num2str(sum(mmm & m100)) ')'], ['JF (' num2str(sum(~mmm & m100)) ')'] },'Interpreter','latex','fontsize',18,'box','off' )
xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$',labFont);
ylabelmine('No. of Hosts',labFont);
set(gca,'yscale','log','fontsize',axFont)
titlemine('TNG100',18)
%text(14.3,400,'TNG100','Fontsize',16,'interpreter','latex')
xlim([10 15.])
ylim([0.5 2000 ])
fname='cjf_jfhost_hist_TNG100';
if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end

myFigure;
histogram(log10(jfStats.M200c(mmm & m50)),linspace(10.4,15,50),'facecolor',cmap(6,:));
hold on;
histogram(log10(jfStats.M200c(~mmm & m50)),linspace(10.4,15,50),'facecolor',cmap(2,:));

histogram(log10(jfStats.M200c(mmm)),linspace(10.4,15,50),'displayStyle','stairs',...
    'edgecolor',cmap(6,:),'linewidth',1.5);
histogram(log10(jfStats.M200c(~mmm)),linspace(10.4,15,50),'displayStyle','stairs',...
    'edgecolor',cmap(2,:),'linewidth',1.5);
%text(10^14.2,1000,num2str(sum(mmm)),'Interpreter','latex','Color',cmap(2,:));
%text(10^14.2,700,num2str(sum(~mmm)),'Interpreter','latex','Color',cmap(1,:));
legend({['Non-JF (' num2str(sum(mmm & m50)) ')'], ['JF (' num2str(sum(~mmm & m50)) ')'] },'Interpreter','latex','fontsize',18,'box','off'   )
xlabelmine('log Host Mass $(M_\mathrm{200,c})\,[\mathrm{M_\odot}]$',labFont);
ylabelmine('No. of Hosts',labFont);
set(gca,'yscale','log','fontsize',axFont)
%text(14.3,400,'TNG50','Fontsize',16,'interpreter','latex')
titlemine('TNG50',18)
xlim([10 15.])
ylim([0.5 2000 ])
fname='cjf_jfhost_hist_TNG50';
if printFlag; printout_fig(gcf,fname,'nopdf','v','printoutdir',outdir); end
