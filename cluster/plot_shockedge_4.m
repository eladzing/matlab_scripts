% plot shockedge array
%  load data
load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat')
%pdir='C:\Users\eladzing\Documents\cluster\printout';
% switch syst
%     case 'win'
%         
%        % load('C:\Users\owner\Ubuntu One\cluster\matlab/mat_files/shockedge.mat');
%         pdir='C:/Users/owner/Ubuntu One/cluster/printout';
%     case 'lin'
%         load('/home/eladzing/Ubuntu One/cluster/matlab/mat_files/shockedge.mat');
%     otherwise
%         error('what is wrong with you? define the system already!');
% end

%printoutdir=sprintf('%s/%s',pdir,'cledge_paper');
%
%define arrays
ll=16;
rv_a1=zeros(ll,1);
mv_a1=zeros(ll,1);
edg1_a1=zeros(ll,2); % by max gradient
edg2_a1=zeros(ll,2); % by maximal entropy
edge_a1=zeros(ll,4)-1; %
edgerange_a1=edg1_a1;

rv_a06=zeros(ll,1);
mv_a06=zeros(ll,1);
edg1_a06=zeros(ll,2); % by max gradient
edg2_a06=zeros(ll,2); % by maximal entropy
edge_a06=zeros(ll,4)-1;


edgerange_a06=edg1_a06;
edgemid_a1=zeros(ll,1);
edgemid_a06=zeros(ll,1);
edgemean_a1=zeros(ll,1);
edgemean_a06=zeros(ll,1);

rv11=[];
mv11=[];
yy11=[];
rv61=[];
mv61=[];
yy61=[];

rv12=[];
mv12=[];
yy12=[];
rv62=[];
mv62=[];
yy62=[];

rv1=[];
mv1=[];
yy1=[];
rv6=[];
mv6=[];
yy6=[];


flg1=true;
flg2=true;
% pre-processing
for i=1:size(shockedge_a1,1)
    name=shockedge_a1{i,1};
    rv_a1(i)=shockedge_a1{i,2};
    mv_a1(i)=shockedge_a1{i,3};
    e1_a1=shockedge_a1{i,4};%./rv_a1(i);
    e2_a1=shockedge_a1{i,5};%./rv_a1(i);
    
    e1_a06=shockedge_a06{i,4};%./rv_a1(i);
    e2_a06=shockedge_a06{i,5};%./rv_a1(i);

    %     edge_a1(i,=cat(2,e1(shockedgeMask1_a1(i,:)),e2(shockedgeMask2_1(i,:)));
%     edg1_a1(i,:)=e1(shockedgeMask1_a1)
%     edg2_a1(i,:)=shockedge_a1{i,5};%./rv_a1(i);
     
        
    rv_a06(i)=shockedge_a06{i,2};
    mv_a06(i)=shockedge_a06{i,3};
      
    ed11=e1_a1(shockedgeMask1_a1(i,:));
    ed061=e1_a06(shockedgeMask1_a06(i,:));
    
    ed12=e2_a1(shockedgeMask2_a1(i,:));
    ed062=e2_a06(shockedgeMask2_a06(i,:));
   
    ed1=cat(2,ed11,ed12);
    ed06=cat(2,ed061,ed062);
   
    edgemean_a1(i)=mean(ed1);
    edgemean_a06(i)=mean(ed06);
    edgerange_a1(i,:)=[min(ed1) max(ed1)];
    edgerange_a06(i,:)=[min(ed06) max(ed06)];
    edgemid_a1(i)=0.5*(min(ed1)+max(ed1));
    edgemid_a06(i)=0.5*(min(ed06)+max(ed06));
    
    for j=1:length(ed11)
        rv11(end+1)=rv_a1(i);
        mv11(end+1)=mv_a1(i);%./1e14;
        yy11(end+1)=ed11(j);
    end
    for j=1:length(ed061)
        rv61(end+1)=rv_a06(i);
        mv61(end+1)=mv_a06(i);%./1e14;
        yy61(end+1)=ed061(j);
    end
    
     for j=1:length(ed12)
        rv12(end+1)=rv_a1(i);
        mv12(end+1)=mv_a1(i);%./1e14;
        yy12(end+1)=ed12(j);
    end
    for j=1:length(ed062)
        rv62(end+1)=rv_a06(i);
        mv62(end+1)=mv_a06(i);%./1e14;
        yy62(end+1)=ed062(j);
    end
    
     for j=1:length(ed1)
        rv1(end+1)=rv_a1(i);
        mv1(end+1)=mv_a1(i);%./1e14;
        yy1(end+1)=ed1(j);
    end
    for j=1:length(ed06)
        rv6(end+1)=rv_a06(i);
        mv6(end+1)=mv_a06(i);%./1e14;
        yy6(end+1)=ed06(j);
    end
end

zz11=yy11./rv11;
zz12=yy12./rv12;
zz61=yy61./rv61;
zz62=yy62./rv62;
zz1=yy1./rv1;
zz6=yy6./rv6;

%% make figures
figure
%hold on 

relist=['u' 'u' 'u' 'r' 'u' 'u' 'u' 'r' 'r' 'u' 'r' 'u' 'r' 'u' 'r' 'u'];
for i=1:16
    switch relist(i)
        case 'u' 
           semilogx(mv_a1(i).*[1 1],edgerange_a1(i,:)./rv_a1(i),'b','linewidth',3.5)
        case 'r'
           semilogx(mv_a1(i).*[1 1],edgerange_a1(i,:)./rv_a1(i),'b','linewidth',3.5)
    end
           if i==1
        hold on
    end
    semilogx(mv_a06(i).*[1 1],edgerange_a06(i,:)./rv_a06(i),'r','linewidth',3.5)
    semilogx([mv_a06(i) mv_a1(i)],[edgemid_a06(i)./rv_a06(i) edgemid_a1(i)./rv_a1(i)],'g--','linewidth',1.5)
end
semilogx(mv11,zz11,'ob','markersize',10.5,'markerfacecolor','auto')
semilogx(mv12,zz12,'db','markersize',10.5,'markerfacecolor','auto')
%plot(mv_a1,edgemean_a1,'--')
semilogx(mv61,zz61,'or','markersize',10.5,'markerfacecolor','auto')
semilogx(mv62,zz62,'dr','markersize',10.5,'markerfacecolor','auto')

hold off
xlim([2e13 3e15])

hl=legend('$z=0$','$z=0.6$');
set(hl,'interpreter','latex','fontsize',12)
xlabelmine('$M_{\mathrm{vir}}\,[\mathrm{M_\odot}]$')
ylabelmine('$R_{edge}/R_{\mathrm{vir}}$')
set(gca,'fontsize',14)
%titlemine('Edge Estimation')
name=sprintf('%s/cledge_vs_mv_all',printoutdir);


%% 
figure
%hold on 

relist=['u' 'u' 'u' 'r' 'u' 'u' 'u' 'r' 'r' 'u' 'r' 'u' 'r' 'u' 'r' 'u'];
for i=1:16
    switch relist(i)
        case 'u' 
           h(1)=semilogx(mv_a1(i).*[1 1],edgerange_a1(i,:)./rv_a1(i),'b','linewidth',2,'DisplayName','$z=0$');
        case 'r'
           semilogx(mv_a1(i).*[1 1],edgerange_a1(i,:)./rv_a1(i),'b','linewidth',2)
    end
    if i==1
        hold on
    end
    h(2)=semilogx(mv_a06(i).*[1 1],edgerange_a06(i,:)./rv_a06(i),'r','linewidth',2,'DisplayName','$z=0.6$');
    semilogx([mv_a06(i) mv_a1(i)],[edgemid_a06(i)./rv_a06(i) edgemid_a1(i)./rv_a1(i)],'--','color',[0 0.75 0],'linewidth',1.5)
end
semilogx(mv11,zz11,'ob','markersize',9,'markerfacecolor','auto','linewidth',1.5)
semilogx(mv12,zz12,'db','markersize',9,'markerfacecolor','auto','linewidth',1.5)
%plot(mv_a1,edgemean_a1,'--')
semilogx(mv61,zz61,'or','markersize',9,'markerfacecolor','auto','linewidth',1.5)
semilogx(mv62,zz62,'dr','markersize',9,'markerfacecolor','auto','linewidth',1.5)

hold off
xlim([3e13 3e15])

hl=legend(h);
set(hl,'interpreter','latex','fontsize',14)
xlabelmine('$M_{vir}\,[\mathrm{M_\odot}]$')
ylabelmine('$R_{edge}/R_{vir}$')
set(gca,'fontsize',14)

set(gcf,'units','centimeters')
pos=get(gcf,'Position');
pos(3)=30;
set(gcf,'Position',pos);
%titlemine('Edge Estimation')
%name=sprintf('%s/cledge_vs_mv_all',printoutdir);
%printout_fig(gcf,'cledge_vs_mv_all')
