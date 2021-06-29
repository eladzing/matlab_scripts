%penetration_depth
clust=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

load('relaxedness.mat')
maxpen=100.*min(pen.penArr(:,:,1),[],2);
maxpen06=100.*min(pen.penArr(:,:,2),[],2);

%% corrections 
maxpen(8)=44; % CL3 
maxpen(13)=41; % cl10
maxpen(15)=44; % cl14
maxpen(16)=6; % cl14

maxpen06(5)=25.5; %cl105 
maxpen06(7)=17; %c107
maxpen06(11)=50; %cl7


rv1=ones(size(clust));
rv06=ones(size(clust));

dx=ones(size(clust));
dy=dx;

dx([10 11 16])=-2;
dx([ 5 7] )=-3;

%dy(7)=-1;

dx=2.0.*dx;
dy=2.0.*dy;

deX=0.05;
deY=0.05./(1+0.5974);

figure('position',[885   244   659   503])

for i=1:length(clust)
    
    new_env(clust(i))
    rv1(i)=get_rvir;
    new_env(clust(i),'a06')
    rv06(i)=get_rvir;
    
    ol='r';
    fc='r';
    disp1='\mathrm{NR}_{z=0}';
    disp2='\mathrm{NR}_{z=0.6}';
    if rlx(i)
        fc='b';
        disp1='\mathrm{R}_{z=0}';
    end
    if rlx06(i)
        ol='b';
        disp2='\mathrm{R}_{z=0.6}';
    end
    
    disp=sprintf('$ %s \\, \\& \\, %s $',disp1,disp2);
    he1=errbar(maxpen(i),maxpen06(i),100.*deY./rv06(i),100.*deY./rv06(i),'-k');
    %errorbar(maxpen(i),maxpen06(i),100.*deY./rv06(i),100.*deY./rv06(i),'.k');
    hold on
    he2=errbar(maxpen(i),maxpen06(i),100.*deX./rv1(i),100.*deX./rv1(i),'-k','horiz');
    set(he1,'linewidth',2)
    set(he2,'linewidth',2)
    %herrorbar(maxpen(i),maxpen06(i),100.*deX./rv1(i),100.*deX./rv1(i),'.k');
    
    
    h(i)=plot(maxpen(i),maxpen06(i),'o','markersize',12','linewidth',3,...
        'markerEdgeColor',ol,'markerFaceColor',fc,'DisplayName',disp);
    %hold on
    
    
    text(maxpen(i)+dx(i),maxpen06(i)+dy(i),num2str(clust(i)),'Fontsize',14)
    
end
hl=legend(h([4 1 8 2]));
set(hl,'Interpreter','latex','Fontsize',14,'Location','NorthWest')
plot([25 25],[0 70],'--k','linewidth',1.5)
plot([0 70],[25 25],'--k','linewidth',1.5)
xlim([0 60])
ylim([0 60])
set(gca,'fontsize',14)
xlabelmine('Max. Penetration at $z=0\, [\%]$',16)
ylabelmine('Max. Penetration at $z=0.6\, [\%]$',16)

printout_fig(gcf,'penetration_relaxedness')