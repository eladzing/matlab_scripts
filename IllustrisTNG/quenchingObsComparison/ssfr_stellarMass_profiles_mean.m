%% plot median ssfr for the three main simulations

%% add in the obs data
global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '/ssfr_rpos_dataGrab.mat'])

global DEFAULT_PRINTOUT_DIR
outDir=[DEFAULT_PRINTOUT_DIR '/obsComp'];

cc=brewermap(8,'Set1');
hf=myFigure;
h=[];
for k=1%:3
    switch k
        case 1
            lt='-';
        case 2
            lt='--';
        case 3
            lt=':';
    end
    
    
    for i=1:4
        
        xx=TNG300(1).starMass(i).xMean;
        y=TNG300(k).starMass(i).yMean;
        y1=log10(y-0.5*TNG300(k).starMass(i).yStanDev);
        y2=log10(y+0.5*TNG300(k).starMass(i).yStanDev);
        yy=log10(TNG300(k).starMass(i).yMean);
%         
%         xx=TNG300(k).starMass(i).xMedian;
%         %y=TNG300(k).starMass(i).yMedian;
%         y1=log10(TNG300(k).starMass(i).yQuarts(2,:));
%         y2=log10(TNG300(k).starMass(i).yQuarts(3,:));
%         yy=log10(TNG300(k).starMass(i).yMedian);
       
        xc=0.05;
        yc=log10(TNG300(k).starMassC(i).yMean);
        ycP=log10(TNG300(k).starMassC(i).yStanDev)-yc;
        ycN=yc-log10(TNG300(k).starMassC(i).yStanDev);
        
        switch i
            case 1
                colo=cc(2,:);
                nam='11-12';
            case 2
                colo=cc(5,:);
                nam='12-13';
                
            case 3
                colo=cc(3,:);
                nam='13-14';
                
            case 4
                colo=cc(1,:);
                nam='14-15';
                
        end
        
        if k==1
            h(i)=plot(xx,yy,'color',colo,'DisplayName',nam,"LineStyle",lt,'linewidth',1.5);
            hold on
            %errorbar(xx,yy,y1,y2,'color',colo)
        patch([xx fliplr(xx)],[y1 fliplr(y2)],colo,...
            'facealpha',0.1,'edgecolor','none')
        elseif k==2 && i==1
            
            h(end+1)=plot(xx,yy,'color',colo,'DisplayName','TNG300-2',"LineStyle",lt);
        elseif k==3 && i==1
            h(end+1)=plot(xx,yy,'color',colo,'DisplayName','TNG300-3',"LineStyle",lt);
        else
            plot(xx,yy,'color',colo,'DisplayName',nam,"LineStyle",lt);
            
            
        end
        
        errorbar(xc,yc,ycN,ycP,'o','color',colo)
        
    end
    
    
    
end
set(gca,'Fontsize',14)

hl=legend(h,'numcolumns',3,'location','northEast','fontsize',14);
xlim([-0.06 2])
ylim([8.5 12.5])
xlabelmine('$r/R_\mathrm{200,c}$');
ylabelmine('log stellar mass');
 titlemine('TNG300 Mean Stellar Mass ' );

 %printout_fig(gcf,'mean_stellarMass_radProf_TNG300','dir',outDir)
 %% sfr 
 
 hf=myFigure;
h=[];
for k=1%:3
    switch k
        case 1
            lt='-';
        case 2
            lt='--';
        case 3
            lt=':';
    end
    
    
    for i=1:4
        
         xx=TNG300(1).ssfrAvg(i).xMean;
         y=TNG300(k).ssfrAvg(i).yMean;
         y1=log10(y-0.5*TNG300(k).ssfrAvg(i).yStanDev);
         y2=log10(y+0.5*TNG300(k).ssfrAvg(i).yStanDev);
         yy=log10(TNG300(k).ssfrAvg(i).yMean);
         
%         xx=TNG300(k).ssfrAvg(i).xMedian;
%         y=TNG300(k).ssfrAvg(i).yMedian;
%         y1=log10(TNG300(k).ssfrAvg(i).yQuarts(2,:));
%         y2=log10(TNG300(k).ssfrAvg(i).yQuarts(3,:));
%         yy=log10(TNG300(k).ssfrAvg(i).yMedian);
       
        xc=0.05;
        yc=log10(TNG300(k).ssfrAvgC(i).yMedian);
        ycP=log10(TNG300(k).ssfrAvgC(i).yQuarts(2,:))-yc;
        ycN=yc-log10(TNG300(k).ssfrAvgC(i).yQuarts(2,:));
        
        switch i
            case 1
                colo=cc(2,:);
                nam='11-12';
            case 2
                colo=cc(5,:);
                nam='12-13';
                
            case 3
                colo=cc(3,:);
                nam='13-14';
                
            case 4
                colo=cc(1,:);
                nam='14-15';
                
        end
        
        if k==1
            h(i)=plot(xx,yy,'-','color',colo,'DisplayName',nam,"LineStyle",lt,'linewidth',1.5);
            hold on
            %errorbar(xx,yy,y1,y2,'color',colo)
        patch([xx fliplr(xx)],[y1 fliplr(y2)],colo,...
            'facealpha',0.1,'edgecolor','none')
        elseif k==2 && i==1
            
            h(end+1)=plot(xx,yy,'color',colo,'DisplayName','TNG300-2',"LineStyle",lt);
        elseif k==3 && i==1
            h(end+1)=plot(xx,yy,'color',colo,'DisplayName','TNG300-3',"LineStyle",lt);
        else
            plot(xx,yy,'color',colo,'DisplayName',nam,"LineStyle",lt);
            
            
        end
        
        errorbar(xc,yc,ycN,ycP,'o','color',colo)
        
    end
    
    for i=5:8
    
    switch i
        case 5
            xx=ssfr_rpos_11_main(2,:);
            yy=ssfr_rpos_11_main(1,:);
            pos=ssfr_rpos_11_top(1,:)-ssfr_rpos_11_main(1,:);
            neg=ssfr_rpos_11_main(1,:)-ssfr_rpos_11_bottom(1,:);
            colo=cc(2,:);
            nam='11-12, Obs';
            
            ccc=ssfr_rpos_11_central(1,2);
            cp=ssfr_rpos_11_central(1,3)-ssfr_rpos_11_central(1,2);
            cm=ssfr_rpos_11_central(1,2)-ssfr_rpos_11_central(1,1);
        case 6
            xx=ssfr_rpos_12_main(2,:);
            yy=ssfr_rpos_12_main(1,:);
            pos=ssfr_rpos_12_top(1,:)-ssfr_rpos_12_main(1,:);
            neg=ssfr_rpos_12_main(1,:)-ssfr_rpos_12_bottom(1,:);
            colo=cc(5,:);
            nam='12-13, Obs';
            
            ccc=ssfr_rpos_12_central(1,2);
            cp=ssfr_rpos_12_central(1,3)-ssfr_rpos_12_central(1,2);
            cm=ssfr_rpos_12_central(1,2)-ssfr_rpos_12_central(1,1);
            
        case 7
            xx=ssfr_rpos_13_main(2,:);
            yy=ssfr_rpos_13_main(1,:);
            pos=ssfr_rpos_13_top(1,:)-ssfr_rpos_13_main(1,:);
            neg=ssfr_rpos_13_main(1,:)-ssfr_rpos_13_bottom(1,:);
            colo=cc(3,:);
            nam='13-14, Obs';
            
            ccc=ssfr_rpos_13_central(1,2);
            cp=ssfr_rpos_13_central(1,3)-ssfr_rpos_13_central(1,2);
            cm=ssfr_rpos_13_central(1,2)-ssfr_rpos_13_central(1,1);
            
        case 8
            xx=ssfr_rpos_14_main(2,:);
            yy=ssfr_rpos_14_main(1,:);
            pos=ssfr_rpos_14_top(1,:)-ssfr_rpos_14_main(1,:);
            neg=ssfr_rpos_14_main(1,:)-ssfr_rpos_14_bottom(1,:);
            colo=cc(1,:);
            nam='14-15, Obs';
            
            ccc=ssfr_rpos_14_central(1,2);
            cp=ssfr_rpos_14_central(1,3)-ssfr_rpos_14_central(1,2);
            cm=ssfr_rpos_14_central(1,2)-ssfr_rpos_14_central(1,1);
            
    end
    
    
    
    
    h(i)=errorbar(xx,yy,neg,pos,'x--','color',colo,...
        'DisplayName',nam);
    
    errorbar(0,ccc,cm,cp,'x','color',colo)
end

    
    
end
set(gca,'Fontsize',14)

hl=legend(h,'numcolumns',2,'location','northEast','fontsize',14);
xlim([-0.06 2])
ylim([-14 -8])
xlabelmine('$r/R_\mathrm{200,c}$');
ylabelmine('log sSFR');
titlemine('TNG300 Mean sSFR' );

printout_fig(gcf,'mean_ssfr_radProf_TNG300','dir',outDir)

%printout_fig(gcf,'median_ssfr_radProf')