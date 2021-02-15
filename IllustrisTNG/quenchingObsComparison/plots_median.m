%% plot median ssfr for the three main simulations 

cc=brewermap(8,'Set1');
hf=myFigure; 


for i=1:4
    
    xx=TNG100(1).ssfrAvg(i).xMedian;
    y1=log10(TNG100(1).ssfrAvg(i).yQuarts(2,:))
    y2=log10(TNG100(1).ssfrAvg(i).yQuarts(3,:))
    yy=log10(TNG100(1).ssfrAvg(i).yMean)
    
   
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
    patch([xx fliplr(xx)],[y1 fliplr(y2)],colo,...
        'facealpha',0.2,'edgecolor','none')
    if i==1; hold on; end
    h(i)=plot(xx,yy,'color',colo,'DisplayName',nam,"LineStyle","-");
   
    
end