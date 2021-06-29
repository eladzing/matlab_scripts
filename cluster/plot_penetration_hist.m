%% make penetration histograms 

penetration_depth_a1;

xb=0.0:0.08:1.0;

[histMax xoutM]=hist(maxpen./100,xb);
[histAvg xoutA]=hist(avgpen./100,xb);
% 
% hf=figure; 
% % 
% bh=[];
% bh(1)=bar(xoutM,histMax,'EdgeColor','none','FaceColor','b','linewidth',2,'DispLayName','Maximal');
% hold on
% bh(2)=bar(xoutA,histAvg,'EdgeColor','none','FaceColor','r','linewidth',2,'DispLayName','Average');
% 
% 
% pH1=get(bh(1),'Children'); %For illustration purpose only.
% pH2=get(bh(2),'Children'); %For illustration purpose only.
% %pH = arrayfun(@(x) allchild(x),bh);
% set(pH1,'FaceAlpha',0.5);
% set(pH2,'FaceAlpha',0.5);
% 
% 
% ylim([0 6.6])
% xlim([0 0.8])
% 
% hl=legend(bh);
% set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest')
% 
% set(gca,'Fontsize',14);
% 
% 
% %xlabel('Penetration Depth [R_{vir}]')
% xlabelmine('Penetration Depth $[\mathrm{R_{vir}}]$')
% ylabelmine('$N$')
% 
% printout_fig(hf,'penen_max_vs_avg');

%% take 2 
hf=figure; 
offSet=0.01;
bh=[];
[xM,yM]=stairs(xoutM,histMax);
bh(1)=area(xM-offSet,yM,'FaceColor','b','DispLayName','Maximal');

hold on
[xA,yA]=stairs(xoutA,histAvg);
bh(2)=area(xA,yA,'FaceColor','r','DispLayName','Average');


pH1=get(bh(1),'Children'); 
pH2=get(bh(2),'Children');
set(pH1,'FaceAlpha',0.7);
set(pH2,'FaceAlpha',0.7);


ylim([0 6.6])
xlim([0 0.8])

hl=legend(bh);
set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest','box','on')
set(gca,'Fontsize',14,'box','on');



%xlabel('Penetration Depth [R_{vir}]')
xlabelmine('Penetration Depth $[\mathrm{R_{vir}}]$')
ylabelmine('$N$')

%printout_fig(hf,'penen_max_vs_avg');
%% relax vs. unrelaxed 
% %xb=0.0:0.1:0.7;
% 
% [histRlx xoutR]=hist(maxpen(rlx)./100,xb);
% [histUrlx xoutU]=hist(maxpen(urlx)./100,xb);
% 
% hf=figure; 
% 
% h=[];
% 
% 
% h(1)=bar(xoutR,histRlx,'stacked','EdgeColor','b','FaceColor','b','linewidth',2,'DispLayName','Relaxed');
% hold on
% h(2)=bar(xoutU,histUrlx,'stacked','EdgeColor','r','FaceColor','r','linewidth',2,'DispLayName','Unrelaxed');
% %alpha(0.4)
% 
% ylim([0 6.6])
% xlim([0 0.8])
% 
% hl=legend(h);
% set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest')
% 
% set(gca,'Fontsize',14);
% 
% xlabelmine('Penetration Depth $[\mathrm{R_{vir}}]$')
% ylabelmine('$N$')
% 
% % %printout_fig(hf,'penen_rlx_vs_urlx');
% 
% %% relax vs. unrelaxed take 2 
% %xb=0.0:0.1:0.7;
% 
% [histRlx xoutR]=hist(maxpen./100,xb);
% [histUrlx xoutU]=hist(maxpen(urlx)./100,xb);
% 
% hf=figure; 
% 
% h=[];
% h(1)=bar(xoutR,histRlx,'EdgeColor','b','FaceColor','none','linewidth',2,'DispLayName','Maximal');
% hold on
% h(2)=bar(xoutU,histUrlx,'EdgeColor','none','FaceColor','r','linewidth',2,'DispLayName','Unrelaxed');
% alpha(0.4)
% 
% ylim([0 6.6])
% xlim([0 0.8])
% 
% hl=legend(h);
% set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest')
% 
% set(gca,'Fontsize',14);
% 
% %xlabel('Penetration Depth [R_{vir}]')
% xlabelmine('Penetration Depth $[\mathrm{R_{vir}}]$')
% ylabelmine('$N$')
% 
% %printout_fig(hf,'penen_rlx_vs_urlx');

%% relax vs. unrelaxed take 3 
%xb=0.0:0.1:0.7;

[histRlx xoutR]=hist(maxpen./100,xb);
[histUrlx xoutU]=hist(maxpen(urlx)./100,xb);



bh=[];
hf=figure; 
[xM,yM]=stairs(xoutM,histMax);

hold on
%bh(2)=area(xM,yM,'FaceColor',[0.7,0.7 1],'DispLayName','Relaxed');


[xU,yU]=stairs(xoutU,histUrlx);
[xR,yR]=stairs(xoutR,histRlx);

%yUt(5:6)=3;yUt(11:12)=4;

%yMt=yM;
%yMt(5:6)=1;yMt(11:12)=3;
yR=yM-yU;
bh(2)=area(xM,yM,'FaceColor',[0.85 0.16 0],'DispLayName','Unrelaxed','LineStyle','none');
bh(3)=area(xR,yR,'FaceColor',[0,0.4 1],'DispLayName','Relaxed','LineStyle','none');
bh(1)=plot(xM,yM,'k','linewidth',3,'DispLayName','Maximal');

% pH1=get(bh(1),'Children'); 
% pH2=get(bh(2),'Children');
% set(pH1,'FaceAlpha',1);
% set(pH2,'FaceAlpha',1);

ylim([0 6.6])
xlim([0 0.8])

hl=legend(bh);
set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest')

set(gca,'Fontsize',14,'box','on');

%xlabel('Penetration Depth [R_{vir}]')
xlabelmine('Penetration Depth $[\mathrm{R_{vir}}]$')
ylabelmine('$N$')

%printout_fig(hf,'penen_rlx_vs_urlx');


%x = 0:1:255;
% figure ('name', 'red') ;
% red = (newImg(:,:,1));
% r = red(:)';
% r = cast(r,'double');
% [graph1,graph2] = hist (r,x);
% bar(graph2,graph1, 'FaceColor', 'r','EdgeColor','r')
% alpha(0.3);
% 
% green = (newImg(:,:,2));
% g = green(:)';
% g = cast(g,'double');
% [graph1,graph2] = hist (g,x);
% bar(graph2,graph1, 'FaceColor', 'b','EdgeColor','b')
% 
% hold off;