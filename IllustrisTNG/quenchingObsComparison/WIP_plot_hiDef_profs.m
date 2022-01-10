hf=figure('color','w','position',figPos);
titlemine('Gal');
h=[];
%yl=[-0.1 0.2];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        profSt=hprofs50.byHostStar4(k).Gal.hiProfD(1,j);
        msk=profSt.binCount>0;
        if any(msk)
            
            
            %          neg=abs(profSt.yQuarts(2,msk)-profSt.yMedian(msk));
            %          pos=abs(profSt.yQuarts(3,msk)-profSt.yMedian(msk));
            px=[profSt.xMedian(msk) fliplr(profSt.xMedian(msk))];
            py=[profSt.yQuarts(2,msk) fliplr(profSt.yQuarts(3,msk))];
            
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            if j==1;hold on;end
            h(end+1)=plot(profSt.xMedian(msk),profSt.yMedian(msk),'-',...
                'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
            %          h(end+1)=errorbar(profSt.xMean,profSt.yMean,profSt.yStanDev./2,'-',...
            %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
            %
            
            %         h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassMedD(1,:,j),'-',...
            %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        end
        
        
        profSt=hprofs100.byHostStar4(k).Gal.hiProfD(1,j);
        msk=profSt.binCount>0;
        if any(msk)
            px=[profSt.xMedian(msk) fliplr(profSt.xMedian(msk))];
            py=[profSt.yQuarts(2,msk) fliplr(profSt.yQuarts(3,msk))];
            
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            if j==1;hold on;end
            h(end+1)=plot(profSt.xMedian(msk),profSt.yMedian(msk),'--',...
                'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        end
        
        %          neg=abs(profSt.yQuarts(2,:)-profSt.yMedian);
        %          pos=abs(profSt.yQuarts(3,:)-profSt.yMedian);
        %         h(end+1)=errorbar(profSt.xMedian,profSt.yMedian,neg,pos,'--',...
        %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        % %          h(end+1)=errorbar(profSt.xMean,profSt.yMean,profSt.yStanDev./2,'--',...
        % %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        %
        %         h(end+1)=plot(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).Gal.hiMassMedD(1,:,j),'--',...
        %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h,'fontsize',14,'location','southeast','interpreter','latex');
        %legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
        %     elseif k==2
        %         legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    %
    
    
    xlim(xl);
    %   ylim(yl);
    xl=xlim;
    yl=ylim;
    
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('HI deficiency',16);
    
    
end



%% average
hf=figure('color','w','position',figPos);
titlemine('Gal');
h=[];
%yl=[-0.1 0.2];
xl=[0 2.4];
for k=1:4
    subplot(2,2,k);
    
    h=[];
    for j=1:4
        
        profSt=hprofs50.byHostStar4(k).Gal.hiProfD(1,j);
        msk=profSt.binCount>0;
        if any(msk)
            
            px=[profSt.xMean(msk) fliplr(profSt.xMean(msk))];
            py=[profSt.yMean(msk)-0.5.*profSt.yStanDev(msk), ...
                fliplr(profSt.yMean(msk)+0.5.*profSt.yStanDev(msk))];
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            if j==1;hold on;end
            h(end+1)=plot(profSt.xMean(msk),profSt.yMean(msk),'-',...
                'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
            %          h(end+1)=errorbar(profSt.xMean,profSt.yMean,profSt.yStanDev./2,'-',...
            %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
            %
            
            %         h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassMedD(1,:,j),'-',...
            %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        end
        
        
        profSt=hprofs100.byHostStar4(k).Gal.hiProfD(1,j);
        msk=profSt.binCount>0;
        if any(msk)
            px=[profSt.xMean(msk) fliplr(profSt.xMean(msk))];
            py=[profSt.yMean(msk)-0.5.*profSt.yStanDev(msk), ...
                fliplr(profSt.yMean(msk)+0.5.*profSt.yStanDev(msk))];
            
            
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            
            h(end+1)=plot(profSt.xMean(msk),profSt.yMean(msk),'--',...
                'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        end
        
        %          neg=abs(profSt.yQuarts(2,:)-profSt.yMedian);
        %          pos=abs(profSt.yQuarts(3,:)-profSt.yMedian);
        %         h(end+1)=errorbar(profSt.xMedian,profSt.yMedian,neg,pos,'--',...
        %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        % %          h(end+1)=errorbar(profSt.xMean,profSt.yMean,profSt.yStanDev./2,'--',...
        % %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        %
        %         h(end+1)=plot(hprofs100.byHostStar4(k).xMed(:,j),hprofs100.byHostStar4(k).Gal.hiMassMedD(1,:,j),'--',...
        %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j});
    end
    
    if k==1
        legend(h,'fontsize',14,'location','southeast','interpreter','latex');
        %         legend(h(1:2:8),'fontsize',14,'location','southeast','interpreter','latex');
        %              elseif k==2
        %                  legend(h(end-1:end),{'TNG50' 'TNG100'},'fontsize',14,'location','northeast','interpreter','latex');
    end
    %
    
    
    xlim(xl);
    %   ylim(yl);
    xl=xlim;
    yl=ylim;
    
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14)
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('HI deficiency',16);
    
    
end

