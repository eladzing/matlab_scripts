
if perlimFlag
    
    global DEFAULT_MATFILE_DIR
    
    tags={'9_95','95_10','10_105','105_11','11_115'};
    tag2={'9'   ,'95'   ,'10'    ,'105'   ,'11'};
    hostTags={'1e12' '5e12' '1e13' '5e13' '1e14' '5e14' '1e15'};
    hostTags2={'12.2' '12.8' '13.2' '13.8' '14.2' '14.8' '15.2'};
    
    baseName='qFrac_Ms100_noRPS_%s_host%s_fgs%s.mat';
    % baseName='qFrac_Ms100_%s_host%s_fgh300_fgs%s.mat';
    
    flags=true(5,7);
    flags(3:5,1)=false;
    flags(5,2:3)=false;
    
    betaInd=1;
    
    for i=1:length(hostTags)
        
        for j=1:length(tags)
            
            if flags(j,i)
                
                
                fname=sprintf(baseName,tags{j},hostTags{i},simEmulate);
                load([DEFAULT_MATFILE_DIR '/' fname ])
                
                qfHost(i).(sprintf('qf%s',tag2{j}))=qFrac;
            end
        end
    end
    
    for i=1:length(hostTags)
        
        for j=1:length(tags)
            
            if j==1
                
                qfHost(i).qfMs1=qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass(betaInd).qf;
                qfHost(i).cntMs1=qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass(betaInd).cnt;
                
                qfHost(i).qfMsF=qfHost(i).(sprintf('qf%s',tag2{j})).msFull(betaInd).qf;
                qfHost(i).cntMsF=qfHost(i).(sprintf('qf%s',tag2{j})).msFull(betaInd).cnt;
                
                qfHost(i).qfMs1P=qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass(betaInd).qfProj;
                qfHost(i).cntMs1P=qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass(betaInd).cntProj;
                
                qfHost(i).qfMsFP=qfHost(i).(sprintf('qf%s',tag2{j})).msFull(betaInd).qfProj;
                qfHost(i).cntMsFP=qfHost(i).(sprintf('qf%s',tag2{j})).msFull(betaInd).cntProj;
                
                qfHost(i).qfMrFP =qfHost(i).(sprintf('qf%s',tag2{j})).massRatSFull(betaInd).qfProj;
                qfHost(i).cntMrFP=qfHost(i).(sprintf('qf%s',tag2{j})).massRatSFull(betaInd).cntProj;
                               
                
                
                qfHost(i).bcen=qfHost(i).(sprintf('qf%s',tag2{j})).msFull(betaInd).bCen-0.25;
            else
                
                if ~isempty(qfHost(i).(sprintf('qf%s',tag2{j})))
                    
                    if length(qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass)>1
                      nb=betaInd;
                    else
                      nb=1;
                    end
                    
                    qfHost(i).qfMs1 = qfHost(i).qfMs1 + qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass(nb).qf;
                    qfHost(i).cntMs1= qfHost(i).cntMs1+ qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass(nb).cnt;
                    
                    qfHost(i).qfMsF = qfHost(i).qfMsF + qfHost(i).(sprintf('qf%s',tag2{j})).msFull(nb).qf;
                    qfHost(i).cntMsF= qfHost(i).cntMsF+ qfHost(i).(sprintf('qf%s',tag2{j})).msFull(nb).cnt;
                    
                    qfHost(i).qfMs1P = qfHost(i).qfMs1P + qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass(nb).qfProj;
                    qfHost(i).cntMs1P= qfHost(i).cntMs1P+ qfHost(i).(sprintf('qf%s',tag2{j})).ms1Pass(nb).cntProj;
                    
                    qfHost(i).qfMsFP = qfHost(i).qfMsFP + qfHost(i).(sprintf('qf%s',tag2{j})).msFull(nb).qfProj;
                    qfHost(i).cntMsFP= qfHost(i).cntMsFP+ qfHost(i).(sprintf('qf%s',tag2{j})).msFull(nb).cntProj;
                    
                    qfHost(i).qfMrFP =  qfHost(i).qfMrFP  + qfHost(i).(sprintf('qf%s',tag2{j})).massRatSFull(nb).qfProj;
                    qfHost(i).cntMrFP=  qfHost(i).cntMrFP + qfHost(i).(sprintf('qf%s',tag2{j})).massRatSFull(nb).cntProj;
                    
                end
            end
            
        end
        
        qfHost(i).qff=qfHost(i).qfMsF./qfHost(i).cntMsF;
        qfHost(i).qf1=qfHost(i).qfMs1./qfHost(i).cntMs1;
        
        qfHost(i).qffP=qfHost(i).qfMsFP./qfHost(i).cntMsFP;
        qfHost(i).qf1P=qfHost(i).qfMs1P./qfHost(i).cntMs1P;
        
        qfHost(i).qfMrP=qfHost(i).qfMrFP./qfHost(i).cntMrFP;
        
        
        
        %     % cut by stellar mass
        %
        %     qfh(i,j)=qfh(i,j)+qfHost(i)
        %     nh(i,j)=nh(i,j)+qfHost(i).
        
        
    end
    
    %% prepare TNG values
    %load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/quenchedFractions_2_TNG_z0.mat')
    load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/quenchedFractions_r200M_TNG_snp99.mat')
    
    msTNG=qFracTNG.stellarMassBins(1:end);
    mhTNG=qFracTNG.hostBins(1:end-1);
    rbTNG=qFracTNG.radialBins(1:4);
    qb100=squeeze(sum(qFracTNG.qbinT100(:,:,1:4),3));
    nb100=squeeze(sum(qFracTNG.nbinT100(:,:,1:4),3));
    
    qb300=squeeze(sum(qFracTNG.qbinT300(:,:,1:4),3));
    nb300=squeeze(sum(qFracTNG.nbinT300(:,:,1:4),3));
    
    qf100=qb100./nb100;
    qf300=qb300./nb300;
    
    
    qbP100=squeeze(sum(qFracTNG.qbinProjT100(:,:,1:4),3));
    nbP100=squeeze(sum(qFracTNG.nbinProjT100(:,:,1:4),3));
    
    qbP300=squeeze(sum(qFracTNG.qbinProjT300(:,:,1:4),3));
    nbP300=squeeze(sum(qFracTNG.nbinProjT300(:,:,1:4),3));
    
    qfP100=qbP100./nbP100;
    qfP300=qbP300./nbP300;
    
    load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/quenchedFractionsCentrals_TNG_z0.mat')
    msCTNG=qFracCenTNG.stellarMassBins(1:end);
    mhCTNG=qFracCenTNG.hostBins(1:end-1);
    
       
    qfC100=qFracCenTNG.qbinT100./qFracCenTNG.nbinT100;
    qfC300=qFracCenTNG.qbinT300./qFracCenTNG.nbinT300;
    
    
    
    %% prepare Wetzel Values
    
    load('/home/zinger/workProjects/matlab_scripts/rpsProject/matFiles/wetzelFull.mat')
    
    
end
%%  plot figures

%% vs stellar mass, cut in halos - Full Pass
cmap=brewermap(8,'Set1');
cmap=cmap([1 2 3 4 5 7 8],:);



% model

% xax(1,:)=qfHost(1).bcen;
% xax(2,:)=qfHost(1).bcen+.5;

xax=qfHost(1).bcen;
xax(end+1)=qfHost(1).bcen(end)+0.5;

xl=[9.0 11.7];
dx=diff(xl)*0.01;
for kk=1:3
    
    switch(kk)
        case 1
            kstart=1;
            kend=2;
            wkstart=1;
            wkend=3;
            
            nc=1;
            
            loc='SouthEast';
            titTag='Low Mass Halos';
            ptag='lohHalos';
        case 2
            kstart=3;
            kend=4;
            wkstart=3;
            wkend=5;
            
            nc=2;
            loc='SouthEast';
            titTag='Intermediate Mass Halos';
            ptag='midHalos';
        case 3
            kstart=5;
            kend=7;
            %kend=6;
            wkstart=5;
            wkend=6;
            nc=3;
            loc='SouthEast';
            titTag='High Mass Halos';
            ptag='hiHalos';
    end
    
    figure('color','w');
    set(gcf,'position',[1432 421 1000 750],'Color','w')
    h=[];
    hcnt=0;
    scount=0;
    for k=kstart:kend
        ltag=sprintf('mod. %s',hostTags2{k});
             ind=find(~isnan(qfHost(k).qffP));
%             yax=cat(1,qfHost(k).qffP,qfHost(k).qffP);
% 
%             hcnt=hcnt+1;
%         
%             hh=plot(xax(:,ind),yax(:,ind),'.-','color',cmap(k,:),...
%                 'DisplayName',ltag,'linewidth',1.5);
%             h(hcnt)=hh(1);
            hcnt=hcnt+1;
            yax=[qfHost(k).qffP(ind) qfHost(k).qffP(ind(end))];
            ddx=(scount-1)*dx;
        h(hcnt)=stairs(xax([ind ind(end)+1])+ddx,yax,'-','color',cmap(k,:),...
            'DisplayName',ltag,'linewidth',1.5);
%         h(hcnt)=stairs(qfHost(k).bcen+0.25,qfHost(k).qffP,'o-','color',cmap(k,:),...
%             'DisplayName',ltag,'linewidth',1.5);
%         
        
        if k==kstart; hold on; end
        scount=scount+1;
    end
    
    % TNG
    
    
    tngTag={'12','12.5','13','13.5','14','14.5','15'};
    scount=0;   
    for k=kstart:kend
        ltag=sprintf('TNG300 %s',tngTag{k});
        ind=find(~isnan(qfP300(:,k)));          
        
        hcnt=hcnt+1;
        
        yax=[qfP300(ind,k); qfP300(ind(end),k)];
            ddx=(scount-1)*dx;
        h(hcnt)=stairs(msTNG([ind ; ind(end)+1])+ddx,yax,'--','color',cmap(k,:),...
            'DisplayName',ltag,'linewidth',1.5);
        
        
%         h(hcnt)=stairs(msTNG+0.25,qfP300(:,k),'-.+','color',cmap(k,:),...
%             'DisplayName',ltag,'linewidth',1.5);
        % if k==1; hold on; end
        scount=scount+1;
    end
%     
%     for k=kstart:kend
%         ltag=sprintf('TNG100 %s',tngTag{k});
%         ind=find(~isnan(qfP100(:,k)));          
%         
%         hcnt=hcnt+1;
%         
%         yax=[qfP100(ind,k); qfP100(ind(end),k)];
%             
%         h(hcnt)=stairs(msTNG([ind ; ind(end)+1]),yax,'--','color',cmap(k,:),...
%             'DisplayName',ltag,'linewidth',1.5);
%         
%         
% %         h(hcnt)=stairs(msTNG+0.25,qfP300(:,k),'-.+','color',cmap(k,:),...
% %             'DisplayName',ltag,'linewidth',1.5);
%         % if k==1; hold on; end
%     end
    
    
    % Wetzel
    fld={'mh12' 'mh125' 'mh13' 'mh135' 'mh14' 'mh145' 'centrals'};
    wTag={'$11.9-12.4$',...
        '$12.4-12.8$',...
        '$12.8-13.3$',...
        '$13.3-13.8$',...
        '$13.8-14.3$',...
        '$14.3-14.8$',...
        'Centrals'};
    
    for k=wkstart:wkend
        
        ltag=sprintf('Wetzel %s',wTag{k});
        
        ms=fig3a.(sprintf('ms_qf_%s',fld{k}))(2,:);
        qf=fig3a.(sprintf('ms_qf_%s',fld{k}))(1,:);
        ypos=abs(fig3a.(sprintf('ms_qf_%s_errP',fld{k}))(1,:)-qf);
        yneg=abs(fig3a.(sprintf('ms_qf_%s_errM',fld{k}))(1,:)-qf);
        hcnt=hcnt+1;
        h(hcnt)=errorbar(ms,qf,yneg,ypos,':','color',cmap(k,:),...
            'Displayname',ltag,'linewidth',1.5);
        %if k==1;hold on;end
    end
    
    xlim(xl)
    ylim([0 1])
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',10,'Location',loc,'NumColumns',nc);
    
    
    xlabelmine('log Stellar Mass $[\mathrm{M_\odot}]$');
    ylabelmine('Quenched Fraction');
    
    titlemine(titTag);
    
    fname=sprintf('qf_comparison_ms_noRPS_fgs%s_%s',simEmulate,ptag);
    %fname=sprintf('qf_comparison_ms_fgh_fgs%s_%s',simEmulate,ptag);
    printout_fig(gcf,fname,'v');
    
end
 
% 
% % pause
%  
% 
% %% sum over all hosts 
% 
% qfCH100=squeeze(sum(qFracCenTNG.qbinT100,2)./sum(qFracCenTNG.nbinT100,2));
% qfCH300=squeeze(sum(qFracCenTNG.qbinT300,2)./sum(qFracCenTNG.nbinT300,2));
% 
% qfH100=squeeze(sum(qb100,2)./sum(nb100,2));
% qfH300=squeeze(sum(qb300,2)./sum(nb300,2));
% qfPH100=squeeze(sum(qbP100,2)./sum(nbP100,2));
% qfPH300=squeeze(sum(qbP300,2)./sum(nbP300,2));
% 
% RomC_data
% 
% 
% figure
% 
% h(1)=stairs(msTNG,[qfPH100 ; qfPH100(end)],'color',cc(1,:),'linestyle','--',...
%     'DisplayName','100 sat','linewidth',1.5);
% hold on
% h(2)=stairs(msTNG,[qfPH300 ; qfPH300(end)],'color',cc(1,:),'linestyle','-',...
% 'DisplayName','300 sat','linewidth',1.5);
% 
% h(3)=stairs(msTNG,[qfCH100 ; qfCH100(end)],'color',cc(2,:),'linestyle','--',...
% 'DisplayName','100 cen','linewidth',1.5);
% 
% h(4)=stairs(msTNG,[qfCH300 ; qfCH300(end)],'color',cc(2,:),'linestyle','-',...
% 'DisplayName','300 cen','linewidth',1.5);
% 
% h(5)=plot(fig3a.ms_qf_centrals(2,:),fig3a.ms_qf_centrals(1,:),...
%     '-dk','DisplayName','Wetzel Cen');
% 
% h(6)=plot(fig3a.ms_qf_mh145(2,:),fig3a.ms_qf_mh145(1,:),...
%     '-ok','DisplayName','Wetzel sat 14.5');
% 
% 
% h(7)=plot(fig14.massBinC(1:end-1)+0.5.*diff(fig14.massBinC),fig14.qfC,'x-',...
%     'color',cc(3,:),'DisplayName','RomC');
% 
% h(8)=plot(fig14.massBin25(1:end-1)+0.5.*diff(fig14.massBin25),fig14.qf25,'x-',...
%     'color',cc(4,:),'DisplayName','Rom25');
% 
% grid
% set(gca,'fontsize',14)
% 
% 
% hl=legend(h);
% set(hl,'Interpreter','latex','Fontsize',12,'location','SouthEast')
% 
% xlabelmine('log Stellar Mass');
% ylabelmine('Quenched Fraction');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % by mass ratio 
% 
% figure
% 
% % model
% h=[];
% edj=-7:1:-1;
% for k=1:7
%     ltag=sprintf('mod. %s',hostTags2{k});
%     
%     h(k)=plot(edj(1:end-1)+0.5,qfHost(k).qfMrP,'-o','color',cmap(k,:),...
%         'DisplayName',ltag,'linewidth',1.5);
%     if k==1; hold on; end
% end
% 
% ylim([0 1])
%     grid
%     
%     hl=legend(h);
%     set(hl,'Interpreter','latex','fontsize',10,'Location',loc,'NumColumns',nc);
%     
%     
%     xlabelmine('log stellar Mass / host haloMass');
%     ylabelmine('Quenched Fraction');
%     
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % %% vs stellar mass, cut in halos  - 1 Pass
% % cmap=brewermap(8,'Set1');
% % cmap=cmap([1 2 3 4 5 7 8],:);
% %
% % figure
% %
% % % model
% % h=[];
% % for k=1:7
% %     ltag=sprintf('mod. %s',hostTags2{k});
% %     h(k)=stairs(qfHost(k).bcen+0.25,qfHost(k).qf1,'-o','color',cmap(k,:),...
% %         'DisplayName',ltag,'linewidth',1.5);
% %     if k==1; hold on; end
% % end
% %
% % % TNG
% % h1=[];
% %
% % tngTag={'12','12.7','13','13.7','14','14.7','15'};
% %
% % for k=1:5
% %     ltag=sprintf('TNG100 %s',tngTag{k});
% %
% %     h(end+1)=plot(msTNG+0.25,qf100(:,k),'--x','color',cmap(k,:),...
% %         'DisplayName',ltag,'linewidth',1.5);
% %     %if k==1; hold on; end
% % end
% %
% % for k=1:7
% %     ltag=sprintf('TNG300 %s',tngTag{k});
% %
% %     h(end+1)=plot(msTNG+0.25,qf300(:,k),'-.+','color',cmap(k,:),...
% %         'DisplayName',ltag,'linewidth',1.5);
% %    % if k==1; hold on; end
% % end
% %
% % % Wetzel
% % fld={'mh12' 'mh125' 'mh13' 'mh135' 'mh14' 'mh145' 'centrals'};
% % wTag={'$11.9-12.4$',...
% %     '$12.4-12.8$',...
% %     '$12.8-13.3$',...
% %     '$13.3-13.8$',...
% %     '$13.8-14.3$',...
% %     '$14.3-14.8$',...
% %     'Centrals'};
% %
% % for k=1:length(fld)-1
% %
% %     ltag=sprintf('Wetzel %s',wTag{k});
% %
% %     ms=fig3a.(sprintf('ms_qf_%s',fld{k}))(2,:);
% %     qf=fig3a.(sprintf('ms_qf_%s',fld{k}))(1,:);
% %     ypos=abs(fig3a.(sprintf('ms_qf_%s_errP',fld{k}))(1,:)-qf);
% %     yneg=abs(fig3a.(sprintf('ms_qf_%s_errM',fld{k}))(1,:)-qf);
% %
% %     h(end+1)=errorbar(ms,qf,yneg,ypos,':','color',cmap(k,:),...
% %         'Displayname',ltag,'linewidth',1.5);
% %     %if k==1;hold on;end
% % end
% %
% %  ylim([0 1])
% % grid
% %
% % hl=legend(h);
% % set(hl,'Interpreter','latex','fontsize',10,'Location','northeast','NumColumns',3);
% %
% %
% % xlabelmine('log Stellar Mass $[\mathrm{M_\odot}]$');
% % ylabelmine('Quenched Fraction');
% %
% %
% %
% %
