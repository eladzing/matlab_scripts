global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/wetzel.mat'])

tags={'9_95','95_10','10_105','105_11','11_115'};


load([DEFAULT_MATFILE_DIR '/qFrac_Ms100_9_95.mat'])
qf9=qFrac;

load([DEFAULT_MATFILE_DIR '/qFrac_Ms100_95_10.mat'])
qf95=qFrac;

load([DEFAULT_MATFILE_DIR '/qFrac_Ms100_10_105.mat'])
qf10=qFrac;

load([DEFAULT_MATFILE_DIR '/qFrac_Ms100_105_11.mat'])
qf105=qFrac;

load([DEFAULT_MATFILE_DIR '/qFrac_Ms100_11_115.mat'])
qf11=qFrac;

clear qFrac


% qff9=qf9.rFull.qfProj./qf9.rFull.cntProj;
% qff95=qf95.rFull.qfProj./qf95.rFull.cntProj;
% qff10=qf10.rFull.qfProj./qf10.rFull.cntProj;
% qff105=qf105.rFull.qfProj./qf105.rFull.cntProj;
% qff11=qf11.rFull.qfProj./qf11.rFull.cntProj;


qfMs1=qf9.ms1Pass.qf+qf95.ms1Pass.qf+qf10.ms1Pass.qf+...
    qf105.ms1Pass.qf+qf11.ms1Pass.qf;

cntMs1=qf9.ms1Pass.cnt+qf95.ms1Pass.cnt+qf10.ms1Pass.cnt+...
    qf105.ms1Pass.cnt+qf11.ms1Pass.cnt;

    
qfMsF=qf9.msFull.qf+qf95.msFull.qf+qf10.msFull.qf+...
    qf105.msFull.qf+qf11.msFull.qf;

cntMsF=qf9.msFull.cnt+qf95.msFull.cnt+qf10.msFull.cnt+...
    qf105.msFull.cnt+qf11.msFull.cnt;

db=diff(smBin);
smb=smBin;  %(1:end-1);


bcen=qf9.msFull.bCen-0.25;

qff=qfMsF./cntMsF;
qff(end+1)=qff(end);
bcen(end+1)=bcen(end)+0.5;

j=3;

%smb(end+1)=smb(end)+db(end);
qff100=qf100(:,j);
qff100(end+1)=qff100(end);
qff300=qf300(:,j);
qff300(end+1)=qff300(end);


figure 
h(3)=stairs(bcen,qff,'color',cmap(3,:),'DisplayName','model','linewidth',1.5);
hold on 


str1=sprintf('$TNG100:[10^{%s}\\,-\\,10^{%s}]$',num2str(hmBin(j-1)),num2str(hmBin(j)));
str3=sprintf('$TNG300:[10^{%s}\\,-\\,10^{%s}]$',num2str(hmBin(j-1)),num2str(hmBin(j)));
    %errorbar(smb,qbin(:,j),db/2,'o')
    
h(1)=stairs(smb,qff100,...
        'color',cmap(1,:),'DisplayName',str1,'linewidth',1.5);
h(2)=stairs(smb,qff300,...
        'color',cmap(2,:),'DisplayName',str3,'linewidth',1.5);
h(4)=plot(Wetzel12_fig14(2,:),Wetzel12_fig14(1,:),'sk','markerfacecolor','k',...
    'DisplayName','Wetzel+12');


ylim([0.5 1])
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',12,'Location','southwest','NumColumns',2);
    
    
    xlabelmine('log Stellar Mass $[\mathrm{M_\odot}]$');
    ylabelmine('Quenched Fraction');
    












%xneg=0.25.*ones(size(bcen));
%xpos=0.25.*ones(size(bcen));

% h(1)=errorbar(bcen,qfMsF./cntMsF,xneg,xpos,'horizontal',...
%     'color',cmap(3,:),'DisplayName','model','linestyle','none','linewidth',1.5);
% hold on 
% 
% j=3;
% str1=sprintf('$100:[10^{%s}\\,-\\,10^{%s}]$',num2str(hmBin(j-1)),num2str(hmBin(j)));
% str3=sprintf('$300:[10^{%s}\\,-\\,10^{%s}]$',num2str(hmBin(j-1)),num2str(hmBin(j)));
%     %errorbar(smb,qbin(:,j),db/2,'o')
%     
% h(2)=errorbar(smb,qf100(:,j),err100(:,j),err100(:,j),db,db,...
%         'color',cmap(1,:),'DisplayName',str1,'linewidth',1.5,'linestyle','none');
% h(3)=errorbar(smb,qf300(:,j),err300(:,j),err300(:,j),db,db,...
%         'color',cmap(2,:),'DisplayName',str3,'linewidth',1.5,'linestyle','none');
% 
% 
% ylim([0.5 1])
%     grid
%     
%     hl=legend(h);
%     set(hl,'Interpreter','latex','fontsize',12,'Location','southwest');%,'NumColumns',3);
%     
%     
%     xlabelmine('log Stellar Mass $[\mathrm{M_\odot}]$');
%     ylabelmine('Quenched Fraction');
%     
% h(1)=plot(bcen,qff9,'-o','DisplayName','9-9.5');
% hold on
% h(2)=plot(bcen,qff95,'-o','DisplayName','9.5-10');
% h(3)=plot(bcen,qff10,'-o','DisplayName','10-10.5');
% h(4)=plot(bcen,qff105,'-o','DisplayName','10.5-11');
% h(5)=plot(bcen,qff11,'-o','DisplayName','11-11.5');

 