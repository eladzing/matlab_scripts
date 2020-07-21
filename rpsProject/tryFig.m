   global DEFAULT_MATFILE_DIR
  
   
    
    load('/home/zinger/workProjects/matlab_scripts/rpsProject/matFiles/wetzelFull.mat')
    
    
   
   tngTag={'12','12.5','13','13.5','14','14.5','15'};

    tags={'9_95','95_10','10_105','105_11','11_115'};
    tag2={'9'   ,'95'   ,'10'    ,'105'   ,'11'};
    hostTags={'1e12' '5e12' '1e13' '5e13' '1e14' '5e14' '1e15'};
    hostTags2={'12.2' '12.8' '13.2' '13.8' '14.2' '14.8' '15.2'};
     fld={'mh12' 'mh125' 'mh13' 'mh135' 'mh14' 'mh145' 'centrals'};
wTag={'$11.9-12.4$',...
        '$12.4-12.8$',...
        '$12.8-13.3$',...
        '$13.3-13.8$',...
        '$13.8-14.3$',...
        '$14.3-14.8$',...
        'Centrals'};
    
    flags=true(5,7);
    flags(3:5,1)=false;
    flags(5,2:3)=false;
    
    betaInd=1;
    
    
    for kk=1:3
        
        switch kk
            case 1 
        baseName='qFrac_Ms100_%s_host%s_fgh300_fgsTNG100.mat';
            case 2
        baseName='qFrac_Ms100_noRPS_%s_host%s_fgsTNG100.mat';
            case 3 
        baseName='qFrac_Ms100_%s_host%s_fgh300_fgs02.mat';

        end
    
       qfHost=struct(); 
    for i=6  %1:length(hostTags)
        
        for j=1:length(tags)
            
            if flags(j,i)
                
                
                fname=sprintf(baseName,tags{j},hostTags{i});%,simEmulate);
                load([DEFAULT_MATFILE_DIR '/' fname ])
                
                qfHost(i).(sprintf('qf%s',tag2{j}))=qFrac;
            end
        end
    end
    
    for i=6  %1:length(hostTags)
        
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
        
        
    end
    
                qfS(kk).qfHost=qfHost;
        
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
    
    
    %% plot 
    
    figure('color','w');
    set(gcf,'position',[1432 421 1000 750],'Color','w')

    h=[];
    hcnt=0;
    scount=0;
    
xl=[9.0 11.7];
dx=diff(xl)*0.01;

xax=qfHost(6).bcen;
xax(end+1)=qfHost(6).bcen(end)+0.5;
cmap=brewermap(8,'Set1');
cmap=cmap([1 2 3 4 5 7 8],:);
    k=6;
    
  
        ltag=sprintf('mod. %s',hostTags2{k});
        ddx=[-1 0 1]*dx;
        
        for kk=1:3
             ind=find(~isnan(qfS(kk).qfHost(k).qffP));
            %hcnt=hcnt+1;
            yax(kk,:)=[qfS(kk).qfHost(k).qffP(ind) qfS(kk).qfHost(k).qffP(ind(end))];
        end    
            
        h(1)=stairs(xax([ind ind(end)+1])+ddx(1),yax(1,:),'-','color',cmap(1,:),...
            'DisplayName','Full Model','linewidth',1.5);

        hold on 

        h(2)=stairs(xax([ind ind(end)+1])+ddx(2),yax(2,:),'-','color',cmap(2,:),...
            'DisplayName','no RPS','linewidth',1.5);

        h(3)=stairs(xax([ind ind(end)+1])+ddx(3),yax(3,:),'-','color',cmap(3,:),...
            'DisplayName','$f_\mathrm{gs}=0.2$','linewidth',1.5);
        
       %% wetzel 
        ltag=sprintf('Wetzel %s',wTag{k});
        
        ms=fig3a.(sprintf('ms_qf_%s',fld{k}))(2,:);
        qf=fig3a.(sprintf('ms_qf_%s',fld{k}))(1,:);
        ypos=abs(fig3a.(sprintf('ms_qf_%s_errP',fld{k}))(1,:)-qf);
        yneg=abs(fig3a.(sprintf('ms_qf_%s_errM',fld{k}))(1,:)-qf);
        hcnt=hcnt+1;
        h(4)=errorbar(ms,qf,yneg,ypos,':','color',cmap(k,:),...
            'Displayname',ltag,'linewidth',1.5);
        
        
        
        
        %% tng 
    
         ltag=sprintf('TNG300 %s',tngTag{k});
        ind=find(~isnan(qfP300(:,k)));          
        
        
        yaxx=[qfP300(ind,k); qfP300(ind(end),k)];
            ddx=(scount-1)*dx;
        h(5)=stairs(msTNG([ind ; ind(end)+1])+ddx,yaxx,'--','color',cmap(k-2,:),...
            'DisplayName',ltag,'linewidth',1.5);
        
        
    ylim([0 1])
    xlim(xl)
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',14,'Location','SouthEast','NumColumns',2);
    
    
    xlabelmine('log Stellar Mass $[\mathrm{M_\odot}]$',16);
    ylabelmine('Quenched Fraction',16);
    set(gca,'fontsize',14)
    
    