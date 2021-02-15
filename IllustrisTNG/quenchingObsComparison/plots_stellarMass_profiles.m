%% plot median ssfr for the three main simulations

cc=brewermap(8,'Set1');
hf=myFigure;
h=[];
for k=1:3
    switch k
        case 1
            lt='-';
        case 2
            lt='--';
        case 3
            lt=':';
    end
    
    
    for i=1:4
        
%         xx=TNG100(1).starMass(i).xMean;
%         y=TNG100(k).starMass(i).yMean;
%         y1=log10(y-0.5*TNG100(k).starMass(i).yStanDev);
%         y2=log10(y+0.5*TNG100(k).starMass(i).yStanDev);
%         yy=log10(TNG100(k).starMass(i).yMean);
%         
        xx=TNG300(k).starMass(i).xMedian;
        y=TNG300(k).starMass(i).yMedian;
        y1=log10(y)-log10(TNG300(k).starMass(i).yQuarts(2,:));
        y2=log10(TNG300(k).starMass(i).yQuarts(3,:))-log10(y);
        yy=log10(TNG300(k).starMass(i).yMedian);
        
        
        
        switch i
            case 1
                colo=cc(2,:);
                nam='11-12, TNG';
            case 2
                colo=cc(5,:);
                nam='12-13, TNG';
                
            case 3
                colo=cc(3,:);
                nam='13-14, TNG';
                
            case 4
                colo=cc(1,:);
                nam='14-15, TNG';
                
        end
        
        if k==1
            h(i)=plot(xx,yy,'color',colo,'DisplayName',nam,"LineStyle",lt);
            hold on
        %patch([xx fliplr(xx)],[y1 fliplr(y2)],colo,...
            errorbar(xx,yy,y1,y2,'color',colo)
        %   'facealpha',0.1,'edgecolor','none')
        elseif k==2 && i==1
            
            h(end+1)=plot(xx,yy,'color',colo,'DisplayName','TNG100-2',"LineStyle",lt);
        elseif k==3 && i==1
            h(end+1)=plot(xx,yy,'color',colo,'DisplayName','TNG100-3',"LineStyle",lt);
        else
            plot(xx,yy,'color',colo,'DisplayName',nam,"LineStyle",lt);
            
            
        end
        hold on
        
        
    end
end
set(gca,'Fontsize',14)

hl=legend(h,'numcolumns',3,'location','northwest','fontsize',14);
xlim([-0.06 1.5])
ylim([9 10.5])
xlabelmine('$r/R_\mathrm{200,c}$');
ylabelmine('log stellar mass');
 titlemine('TNG100 Mean Stellar Mass resolution dependence' );

 
 %% sfr 