xx=categorical({'$\ge 1 R_{200}$' '$\ge1.5 R_{200}$' '$\ge 2 R_{200}$'});
xx=reordercats(xx,{'$\ge 1 R_{200}$' '$\ge1.5 R_{200}$' '$\ge 2 R_{200}$'});
%% 
maskJF=objectTable.scoreWeighted'>=0.8;
mask50=objectTable.sim'=="TNG50";
dmask1=galProps.rpos./galProps.hostR200c>=1;
dmask15=galProps.rpos./galProps.hostR200c>=1.5;
dmask2=galProps.rpos./galProps.hostR200c>=2;

%jf=[24.5 13.2 ; 7.2 5.2 ; 3.6 2.8];
%njf=[49 47.2 ; 26.5 28.3 ; 12.4 15.2]; 

jf=[sum(dmask1 & maskJF & mask50)./sum(maskJF & mask50) sum(dmask1 & maskJF & ~mask50)./sum(maskJF & ~mask50); ...
    sum(dmask15 & maskJF & mask50)./sum(maskJF & mask50) sum(dmask15 & maskJF & ~mask50)./sum(maskJF & ~mask50); ...
    sum(dmask2 & maskJF & mask50)./sum(maskJF & mask50) sum(dmask2 & maskJF & ~mask50)./sum(maskJF & ~mask50)];

njf=[sum(dmask1 & ~maskJF & mask50)./sum(~maskJF & mask50) sum(dmask1 & ~maskJF & ~mask50)./sum(~maskJF & ~mask50); ...
    sum(dmask15 & ~maskJF & mask50)./sum(~maskJF & mask50) sum(dmask15 & ~maskJF & ~mask50)./sum(~maskJF & ~mask50); ...
    sum(dmask2 & ~maskJF & mask50)./sum(~maskJF & mask50) sum(dmask2 & ~maskJF & ~mask50)./sum(~maskJF & ~mask50)];

jf=round(jf.*100,1);
njf=round(njf.*100,1);
%% plot
hf=myFigure;

cols=brewermap(12,'Paired');
%% 

w1=0.8;
w2=0.6;
figure(hf)
hold off
b1=bar(xx,njf,w1,'facecolor',cols(6,:));
hold on
b2=bar(xx,jf,w2,'facecolor',cols(2,:));

dx=0.0;
dy=-0.3;
xtips50 = b1(1).XEndPoints-dx;
xtips100 = b1(2).XEndPoints+dx;
ytips50 = b1(1).YEndPoints+dy;
ytips100 = b1(2).YEndPoints+dy;
ytips=max(ytips50,ytips100)+ dy;
 
labels50 = ["TNG50" "TNG50" "TNG50" ];
labels100 = ["TNG100" "TNG100" "TNG100" ];
text(xtips50,ytips50,labels50,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Interpreter','latex','fontsize',16,'rotation',270)
text(xtips100,ytips100,labels100,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Interpreter','latex','fontsize',16,'rotation',270)
b2(1).DisplayName='JF Population Fraction';
b1(1).DisplayName='Non-JF Population Fraction';

legend([b2(1) b1(1)],'Interpreter','latex','fontsize',16,'location','northeast')

%xlim([0.5 3.5])
ylim([0 52])
ylabelmine('Percentage')
set(gca,'TickLabelInterpreter','latex','fontsize',18);
%% print fig

global DEFAULT_PRINTOUT_DIR
outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/paper'];
outname='cjf_distance_popFraction';
if printFlag; printout_fig(gcf,outname,'nopdf','v','printoutdir',outdir); end


