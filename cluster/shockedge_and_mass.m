% plot shockedge array
%  load data

switch syst
    case 'win'
        load('C:\Users\owner\Ubuntu One\cluster\matlab/mat_files/shockedge.mat');
        pdir='C:/Users/owner/Ubuntu One/cluster/printout';
    case 'lin'
        load('/home/eladzing/Ubuntu One/cluster/matlab/mat_files/shockedge.mat');
    otherwise
        error('what is wrong with you? define the system already!');
end
printoutdir=sprintf('%s/%s',pdir,'cledge_paper');
%
%define arrays
ll=16;
rv_a1=zeros(ll,1);
mv_a1=zeros(ll,1);
edg1_a1=zeros(ll,2); % by max gradient
edg2_a1=zeros(ll,2); % by maximal entropy
edge_a1=zeros(ll,4)-1; %
edgerange_a1=edg1_a1;
mmrange_a1=edgerange_a1;

rv_a06=zeros(ll,1);
mv_a06=zeros(ll,1);
edg1_a06=zeros(ll,2); % by max gradient
edg2_a06=zeros(ll,2); % by maximal entropy
edge_a06=zeros(ll,4)-1;


edgerange_a06=edg1_a06;
mmrange_a106=edgerange_a06;
edgemid_a1=zeros(ll,1);
edgemid_a06=zeros(ll,1);
edgemean_a1=zeros(ll,1);
edgemean_a06=zeros(ll,1);

mshk11=zeros(ll,4,2);
mshk12=zeros(ll,4,2);
mshk061=zeros(ll,4,2);
mshk062=zeros(ll,4,2);

rv11=[];
mv11=[];
yy11=[];
mm11=[];

rv61=[];
mv61=[];
yy61=[];
mm61=[];

rv12=[];
mv12=[];
yy12=[];
mm12=[];

rv62=[];
mv62=[];
yy62=[];
mm62=[];



rv1=[];
mv1=[];
yy1=[];
rv6=[];
mv6=[];
yy6=[];

mdum=[];

flg1=true;
flg2=true;
% pre-processing
for i=1:size(shockedge_a1,1)
    name=shockedge_a1{i,1};
    rv_a1(i)=shockedge_a1{i,2};
    mv_a1(i)=shockedge_a1{i,3};
    edg1_a1(i,:)=shockedge_a1{i,4};%./rv_a1(i);
    edg2_a1(i,:)=shockedge_a1{i,5};%./rv_a1(i);
    
    rv_a06(i)=shockedge_a06{i,2};
    mv_a06(i)=shockedge_a06{i,3};
    edg1_a06(i,:)=shockedge_a06{i,4};%./rv_a06(i);
    edg2_a06(i,:)=shockedge_a06{i,5};%./rv_a06(i);
    
    
    
    switch name
        case('CL101')
            if flg1 ; edge_a1(i,1)=edg1_a1(i,1); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1:2)=edg1_a06(i,:); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL102')
            if flg1 ; edge_a1(i,1:2)=edg1_a1(i,:); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3:4)=edg2_a06(i,:); end
            
        case('CL103')
            if flg1 ; edge_a1(i,1:2)=edg1_a1(i,:); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL104')
            if flg1 ; edge_a1(i,1)=edg2_a1(i,2); end
            if flg2 ; edge_a1(i,3)=edg1_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL105')
            if flg1 ; edge_a1(i,1)=edg1_a1(i,1); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL106')
            if flg1 ; edge_a1(i,1)=edg1_a1(i,1); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
            
        case('CL107')
            if flg1 ; edge_a1(i,1:2)=edg1_a1(i,1:2); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1:2)=edg1_a06(i,1:2); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL3')
            if flg1 ; edge_a1(i,1)=edg1_a1(i,1); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL5')
            if flg1 ; edge_a1(i,1)=edg1_a1(i,1); end
            if flg2 ; edge_a1(i,3:4)=edg2_a1(i,1:2); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL6')
            if flg1 ; edge_a1(i,1:2)=edg1_a1(i,1:2); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1:2)=edg1_a06(i,:); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL7')
            if flg1 ; edge_a1(i,1)=edg1_a1(i,1); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL9')
            if flg1 ; edge_a1(i,1)=edg1_a1(i,1); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL10')
            if flg1 ; edge_a1(i,1)=edg1_a1(i,1); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,1); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL11')
            if flg1 ; edge_a1(i,1:2)=edg1_a1(i,:); end
            if flg2 ; edge_a1(i,3:4)=edg2_a1(i,:); end
            
            if flg1 ; edge_a06(i,1)=edg1_a06(i,1); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
            
        case('CL14')
            if flg1 ; edge_a1(i,1:2)=edg1_a1(i,:); end
            if flg2 ; edge_a1(i,3:4)=edg2_a1(i,:); end
            
            if flg1 ; edge_a06(i,1:2)=edg1_a06(i,:); end
            if flg2 ; edge_a06(i,3:4)=edg2_a06(i,:); end
            
        case('CL24')
            if flg1 ; edge_a1(i,1:2)=edg1_a1(i,:); end
            if flg2 ; edge_a1(i,3)=edg2_a1(i,2); end
            
            if flg1 ; edge_a06(i,1:2)=edg1_a06(i,1:2); end
            if flg2 ; edge_a06(i,3)=edg2_a06(i,1); end
    end
    
    ed11=edge_a1(i,1:2);
    ed11=ed11(ed11>0);
    
    ed061=edge_a06(i,1:2);
    ed061=ed061(ed061>0);
    
    ed12=edge_a1(i,3:4);
    ed12=ed12(ed12>0);
    
    ed062=edge_a06(i,3:4);
    ed062=ed062(ed062>0);
    
    ed1=edge_a1(i,:);
    ed1=ed1(ed1>0);
    
    ed06=edge_a06(i,:);
    ed06=ed06(ed06>0);
    
     
%     
     
    
    
    edgemean_a1(i)=mean(ed1);
    edgemean_a06(i)=mean(ed06);
    edgerange_a1(i,:)=[min(ed1) max(ed1)];
    edgerange_a06(i,:)=[min(ed06) max(ed06)];
    edgemid_a1(i)=0.5*(min(ed1)+max(ed1));
    edgemid_a06(i)=0.5*(min(ed06)+max(ed06));
    
    new_env(name,'csf','a1');
    
    [mshk11(i,1,1:length(ed11)),mshk11(i,2,1:length(ed11)),mshk11(i,3,1:length(ed11))]=read_Mass_Profiles(ed11);
    [mshk12(i,1,1:length(ed12)),mshk12(i,2,1:length(ed12)),mshk12(i,3,1:length(ed12))]=read_Mass_Profiles(ed12);
    
    [a b ~]=read_Mass_Profiles(ed11);
    mt11=a+b;   %read_MTOT_Profile(ed11);
    [a b ~]=read_Mass_Profiles(ed12);
    mt12=a+b; %read_MTOT_Profile(ed12);
    
    new_env(name,'csf','a06');
    [mshk061(i,1,1:length(ed061)),mshk061(i,2,1:length(ed061)),mshk061(i,3,1:length(ed061))]=read_Mass_Profiles(ed061);
    [mshk062(i,1,1:length(ed062)),mshk062(i,2,1:length(ed062)),mshk062(i,3,1:length(ed062))]=read_Mass_Profiles(ed062);
    
    [a b ~]=read_Mass_Profiles(ed061);
    mt061=a+b;%read_MTOT_Profile(ed061);
    [a b ~]=read_Mass_Profiles(ed062);
    mt062=a+b;%read_MTOT_Profile(ed062);
     
    mshk11(i,4,:)=sum(mshk11(i,:,:),2);
    mshk12(i,4,:)=sum(mshk12(i,:,:),2);
    mshk061(i,4,:)=sum(mshk061(i,:,:),2);
    mshk062(i,4,:)=sum(mshk062(i,:,:),2);
    
    mm=cat(2,mt11,mt12);
    mmrange_a1(i,:)=[min(mm) max(mm)];
    mm=cat(2,mt061,mt062);
    mmrange_a06(i,:)=[min(mm) max(mm)];
    
    for j=1:length(ed11)
        rv11(end+1)=rv_a1(i);
        mv11(end+1)=mv_a1(i);%./1e14;
        yy11(end+1)=ed11(j);
        mm11(end+1)=mt11(j);
    end
    for j=1:length(ed061)
        rv61(end+1)=rv_a06(i);
        mv61(end+1)=mv_a06(i);%./1e14;
        yy61(end+1)=ed061(j);
        mm61(end+1)=mt061(j);
    end
    
     for j=1:length(ed12)
        rv12(end+1)=rv_a1(i);
        mv12(end+1)=mv_a1(i);%./1e14;
        yy12(end+1)=ed12(j);
        mm12(end+1)=mt12(j);
    end
    for j=1:length(ed062)
        rv62(end+1)=rv_a06(i);
        mv62(end+1)=mv_a06(i);%./1e14;
        yy62(end+1)=ed062(j);
        mm62(end+1)=mt062(j);
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

xx11=mm11./mv11;
xx12=mm12./mv12;
xx61=mm61./mv61;
xx62=mm62./mv62;

zz1=yy1./rv1;
zz6=yy6./rv6;

%% make figures
figure

%hold on 


for i=1:16
    semilogx(mv_a1(i).*[1 1],edgerange_a1(i,:)./rv_a1(i),'b','linewidth',3.5)
    if i==1
        hold on
    end
    semilogx(mv_a06(i).*[1 1],edgerange_a06(i,:)./rv_a06(i),'r','linewidth',3.5)
    
end
semilogx(mv11,zz11,'ob','markersize',10.5,'markerfacecolor','auto')
semilogx(mv12,zz12,'db','markersize',10.5,'markerfacecolor','auto')
%plot(mv_a1,edgemean_a1,'--')
semilogx(mv61,zz61,'or','markersize',10.5,'markerfacecolor','auto')
semilogx(mv62,zz62,'dr','markersize',10.5,'markerfacecolor','auto')

hold off
xlim([6e13 3e15])

hl=legend('$z=0$','$z=0.6$');
set(hl,'interpreter','latex','fontsize',12)
xlabelmine('$M_{vir}\,[\mathrm{M_\odot}]$')
ylabelmine('$R_{edge}/R_{vir}$')
titlemine('Edge Estimation')
%name=sprintf('%s/cledge_vs_mv_all',printoutdir)
%exportfig(gcf,sprintf('%s.png',name),'format','png');
%exportfig(gcf,sprintf('%s.eps',name));
%------------------------
figure
for i=1:16
    semilogx(mv_a1(i).*[1 1],mmrange_a1(i,:)./mv_a1(i),'b','linewidth',3.5)
    if i==1
        hold on
    end
    semilogx(mv_a06(i).*[1 1],mmrange_a06(i,:)./mv_a06(i),'r','linewidth',3.5)
   
end

semilogx(mv11,xx11,'ob','markersize',10.5,'markerfacecolor','auto')
hold on
semilogx(mv12,xx12,'db','markersize',10.5,'markerfacecolor','auto')
%plot(mv_a1,edgemean_a1,'--')
semilogx(mv61,xx61,'or','markersize',10.5,'markerfacecolor','auto')
semilogx(mv62,xx62,'dr','markersize',10.5,'markerfacecolor','auto')

hold off
xlim([6e13 3e15])

hl=legend('$z=0$','$z=0.6$');
set(hl,'interpreter','latex','fontsize',12)
xlabelmine('$M_{vir}\,[\mathrm{M_\odot}]$')
ylabelmine('$M_{edge}/M_{vir}$')
titlemine('Edge Estimation')
%name=sprintf('%s/cledge_vs_mv_all',printoutdir)
%exportfig(gcf,sprintf('%s.png',name),'format','png');
%exportfig(gcf,sprintf('%s.eps',name));
%% ----------------------------------------------------------------



% figure
% errorbar(mv_a1./1e14,edgemean_a1,edgemean_a1-edgerange_a1(:,1),abs(edgemean_a1-edgerange_a1(:,2)),'o','markersize',8)
%
% figure
% errorbar(mv_a06./1e14,edgemean_a06,edgemean_a06-edgerange_a06(:,1),abs(edgemean_a06-edgerange_a06(:,2)),'o','markersize',8)
% xlabelmine('$M_{vir}\,[10^{14}\mathrm{M_\odot}]$')
% ylabelmine('$R_{edge}/R_{vir}$') titlemine('Edge Estimation at $z=0.6$')
% exportfig(gcf,'edge_vs_mv_a06.png','format','png');
% exportfig(gcf,'edge_vs_mv_a06.eps')

% figure
% for i=1:16
%     semilogx(rv_a1(i).*[1 1],edgerange_a1(i,:),'b','linewidth',3.5)
%     if i==2
%         hold on
%     end
%     semilogx(rv_a06(i).*[1 1],edgerange_a06(i,:),'r','linewidth',3.5)
%     
% end
% semilogx(rv1,yy1,'ob','markersize',10.5,'markerfacecolor','auto')
% semilogx(rv6,yy6,'or','markersize',10.5,'markerfacecolor','auto')
% hold off
% hl=legend('$z=0$','$z=0.6$');
% set(hl,'interpreter','latex','fontsize',12)
% xlabelmine('$R_{vir}\,[\mathrm{Mpc}]$')
% ylabelmine('$R_{edge}/R_{vir}$')
% titlemine('Edge Estimation at $z=0.6$')
% name=sprintf('%s/cledge_vs_rv_a06',printoutdir)
% %exportfig(gcf,sprintf('%s.png',name),'format','png');
% %exportfig(gcf,sprintf('%s.eps',name));
% 

% 
% 
% 
% % figure
% % for i=1:16
% %     plot(rv_a06(i).*[1 1],edgerange_a06(i,:),'linewidth',3.5)
% %     hold on
% % end
% % plot(rv6,yy6,'ob','markersize',10.5,'markerfacecolor','auto')
% % hold off
% % xlabelmine('$R_{vir}\,[\mathrm{Mpc}]$')
% %  ylabelmine('$R_{edge}/R_{vir}$')
% %  titlemine('Edge Estimation at $z=0.6$')
% % % exportfig(gcf,'edge_vs_mv_a1.png','format','png');
% % % exportfig(gcf,'edge_vs_mv_a1.eps')
% %
% % figure
% % for i=1:16
% %     plot(mv_a06(i)/1e14.*[1 1],edgerange_a06(i,:),'linewidth',3.5)
% %     hold on
% % end
% % plot(mv6,yy6,'ob','markersize',10.5,'markerfacecolor','auto')
% % hold off
% % xlabelmine('$M_{vir}\,[10^{14}\mathrm{M_\odot}]$')
% %  ylabelmine('$R_{edge}/R_{vir}$')
% %  titlemine('Edge Estimation at $z=0.6$')
% % % exportfig(gcf,'edge_vs_mv_a1.png','format','png');
% % % exportfig(gcf,'edge_vs_mv_a1.eps')
