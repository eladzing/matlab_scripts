% plot shockedge array 
%  load data 
load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat')
% switch syst
%     case 'win'
%         load('C:\Users\owner\Ubuntu One\cluster\matlab/mat_files/shockedge.mat');
%     case 'lin'
%         load('/home/eladzing/Ubuntu One/cluster/matlab/mat_files/shockedge.mat');
%     otherwise
%         error('what is wrong with you? define the system already!');
% end
%
%define arrays
ll=16;
rv_a1=zeros(ll,1);
mv_a1=zeros(ll,1);
edg1_a1=zeros(ll,2); % by max gradient 
edg2_a1=zeros(ll,2); % by maximal entropy 
edge_a1=zeros(ll,4)-1; % by maximal entropy 
edgerange_a1=edg1_a1;

rv_a06=zeros(ll,1);
mv_a06=zeros(ll,1);
edg1_a06=zeros(ll,2); % by max gradient 
edg2_a06=zeros(ll,2); % by maximal entropy 
edge_a06=zeros(ll,4)-1; % by maximal entropy 
edgerange_a06=edg1_a06;
edgemid_a1=zeros(ll,1);
edgemid_a06=zeros(ll,1);
edgemean_a1=zeros(ll,1);
edgemean_a06=zeros(ll,1);

rv1=[];
mv1=[];
yy1=[];
rv6=[];
mv6=[];
yy6=[];

% pre-processing 
for i=1:size(shockedge_a1,1)
    name=shockedge_a1{i,1};
    rv_a1(i)=shockedge_a1{i,2};
    mv_a1(i)=shockedge_a1{i,3};
    edg1_a1(i,:)=shockedge_a1{i,4}./rv_a1(i);
    edg2_a1(i,:)=shockedge_a1{i,5}./rv_a1(i);
    
    rv_a06(i)=shockedge_a06{i,2};
    mv_a06(i)=shockedge_a06{i,3};
    edg1_a06(i,:)=shockedge_a06{i,4}./rv_a06(i);
    edg2_a06(i,:)=shockedge_a06{i,5}./rv_a06(i);
    
    switch name
        case('CL101')
            edge_a1(i,1)=edg1_a1(i,1); 
            edge_a1(i,2)=edg2_a1(i,1);  
            
            edge_a06(i,1:2)=edg1_a06(i,:); 
            edge_a06(i,3:4)=edg2_a06(i,:);  
            
        case('CL102')
            edge_a1(i,1:2)=edg1_a1(i,:); 
            edge_a1(i,3)=edg2_a1(i,1);  
            
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2:3)=edg2_a06(i,:);  
        
        case('CL103')
            edge_a1(i,1:2)=edg1_a1(i,:); 
            edge_a1(i,3)=edg2_a1(i,1); 
            
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        case('CL104')
            edge_a1(i,1)=edg2_a1(i,2); 
            edge_a1(i,2)=edg1_a1(i,1);  
        
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        case('CL105')
            edge_a1(i,1)=edg1_a1(i,1); 
            edge_a1(i,2)=edg2_a1(i,1);  
            
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        case('CL106')
            edge_a1(i,1)=edg1_a1(i,1); 
            edge_a1(i,2)=edg2_a1(i,1);  
            
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        
        case('CL107')
            edge_a1(i,1:2)=edg1_a1(i,1:2); 
            edge_a1(i,3)=edg2_a1(i,1);  
            
            edge_a06(i,1:2)=edg1_a06(i,1:2); 
            edge_a06(i,3)=edg2_a06(i,1);  
            
        case('CL3')
            edge_a1(i,1)=edg1_a1(i,1); 
            edge_a1(i,2)=edg2_a1(i,1);  
        
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1); 
            
        case('CL5')
            edge_a1(i,1)=edg1_a1(i,1); 
            edge_a1(i,2:3)=edg2_a1(i,1:2);  
        
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        case('CL6')
            edge_a1(i,1:2)=edg1_a1(i,1:2); 
            edge_a1(i,3)=edg2_a1(i,1);  
        
            edge_a06(i,1:2)=edg1_a06(i,:); 
            edge_a06(i,3)=edg2_a06(i,1);  
            
        case('CL7')
            edge_a1(i,1)=edg1_a1(i,1); 
            edge_a1(i,2)=edg2_a1(i,1);  
        
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        case('CL9')
            edge_a1(i,1)=edg1_a1(i,1); 
            edge_a1(i,2)=edg2_a1(i,1);  
            
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        case('CL10')
            edge_a1(i,1)=edg1_a1(i,1); 
            edge_a1(i,2)=edg2_a1(i,1);  
        
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        case('CL11')
            edge_a1(i,1:2)=edg1_a1(i,:); 
            edge_a1(i,3:4)=edg2_a1(i,:);  
        
            edge_a06(i,1)=edg1_a06(i,1); 
            edge_a06(i,2)=edg2_a06(i,1);  
            
        case('CL14')
            edge_a1(i,1:2)=edg1_a1(i,:); 
            edge_a1(i,3:4)=edg2_a1(i,:);
        
            edge_a06(i,1:2)=edg1_a06(i,:); 
            edge_a06(i,3:4)=edg2_a06(i,:);  
            
        case('CL24')
            edge_a1(i,1:2)=edg1_a1(i,:); 
            edge_a1(i,3)=edg2_a1(i,2);
            
            edge_a06(i,1:2)=edg1_a06(i,1:2); 
            edge_a06(i,3)=edg2_a06(i,1);  
    end
    ed1=edge_a1(i,:);
    ed1=ed1(ed1>0);
    
    ed06=edge_a06(i,:);
    ed06=ed06(ed06>0);
    
    edgemean_a1(i)=mean(ed1);
    edgemean_a06(i)=mean(ed06);
    edgerange_a1(i,:)=[min(ed1) max(ed1)];
    edgerange_a06(i,:)=[min(ed06) max(ed06)];
    edgemid_a1(i)=0.5*(min(ed1)+max(ed1));
    edgemid_a06(i)=0.5*(min(ed06)+max(ed06));

    for j=1:length(ed1)
        rv1(end+1)=rv_a1(i);
        mv1(end+1)=mv_a1(i)./1e14;
        yy1(end+1)=ed1(j);
    end
    for j=1:length(ed06)
        rv6(end+1)=rv_a06(i);
        mv6(end+1)=mv_a06(i)./1e14;
        yy6(end+1)=ed06(j);
    end

end   
 
% make figures

% figure 
% errorbar(mv_a1./1e14,edgemean_a1,edgemean_a1-edgerange_a1(:,1),abs(edgemean_a1-edgerange_a1(:,2)),'o','markersize',8)
% 
% figure
% errorbar(mv_a06./1e14,edgemean_a06,edgemean_a06-edgerange_a06(:,1),abs(edgemean_a06-edgerange_a06(:,2)),'o','markersize',8)
% xlabelmine('$M_{vir}\,[10^{14}\mathrm{M_\odot}]$')
% ylabelmine('$R_{edge}/R_{vir}$') titlemine('Edge Estimation at $z=0.6$')
% exportfig(gcf,'edge_vs_mv_a06.png','format','png');
% exportfig(gcf,'edge_vs_mv_a06.eps')

figure
for i=1:16
    plot(rv_a1(i).*[1 1],edgerange_a1(i,:),'b','linewidth',3.5) 
    hold on
    plot(rv_a06(i).*[1 1],edgerange_a06(i,:),'r','linewidth',3.5)
    
end
plot(rv1,yy1,'ob','markersize',10.5,'markerfacecolor','auto')
plot(rv6,yy6,'or','markersize',10.5,'markerfacecolor','auto')
hold off
xlabelmine('$R_{vir}\,[\mathrm{Mpc}]$')
 ylabelmine('$R_{edge}/R_{vir}$')
 titlemine('Edge Estimation at $z=0$')
% exportfig(gcf,'edge_vs_mv_a1.png','format','png');
% exportfig(gcf,'edge_vs_mv_a1.eps')
 
figure
for i=1:16
    plot(mv_a1(i)/1e14.*[1 1],edgerange_a1(i,:),'b','linewidth',3.5)
    hold on
    plot(mv_a06(i)/1e14.*[1 1],edgerange_a06(i,:),'r','linewidth',3.5)
    
end
plot(mv1,yy1,'ob','markersize',10.5,'markerfacecolor','auto')
plot(mv6,yy6,'or','markersize',10.5,'markerfacecolor','auto')
hold off
xlabelmine('$M_{vir}\,[10^{14}\mathrm{M_\odot}]$')
 ylabelmine('$R_{edge}/R_{vir}$')
 titlemine('Edge Estimation at $z=0$')
% exportfig(gcf,'edge_vs_mv_a1.png','format','png');
% exportfig(gcf,'edge_vs_mv_a1.eps')

% figure
% for i=1:16
%     plot(rv_a06(i).*[1 1],edgerange_a06(i,:),'linewidth',3.5)
%     hold on
% end
% plot(rv6,yy6,'ob','markersize',10.5,'markerfacecolor','auto')
% hold off
% xlabelmine('$R_{vir}\,[\mathrm{Mpc}]$')
%  ylabelmine('$R_{edge}/R_{vir}$')
%  titlemine('Edge Estimation at $z=0.6$')
% % exportfig(gcf,'edge_vs_mv_a1.png','format','png');
% % exportfig(gcf,'edge_vs_mv_a1.eps')
%  
% figure
% for i=1:16
%     plot(mv_a06(i)/1e14.*[1 1],edgerange_a06(i,:),'linewidth',3.5)
%     hold on
% end
% plot(mv6,yy6,'ob','markersize',10.5,'markerfacecolor','auto')
% hold off
% xlabelmine('$M_{vir}\,[10^{14}\mathrm{M_\odot}]$')
%  ylabelmine('$R_{edge}/R_{vir}$')
%  titlemine('Edge Estimation at $z=0.6$')
% % exportfig(gcf,'edge_vs_mv_a1.png','format','png');
% % exportfig(gcf,'edge_vs_mv_a1.eps')
