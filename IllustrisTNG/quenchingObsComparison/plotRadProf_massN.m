
function plotRadProf(profStr,varargin)  %,type,sim,model,printFlag)

printFlag=false;
type='none';
sim='none';
model=0;

outDir='C:\Users\eladz\Documents\workProjects\IllustrisTNG\printout\obsComp';

%% parse arguments
i=1;
while(i<=length(varargin))
    
    switch(lower(varargin{i}))
        case{'print','plot','save'}
            printFlag=true;
        case{50,'50','tng50'}
            sim='TNG50';
        case{100,'100','tng100'}
            sim='TNG100';
        case{1,'br'}
            model=1;
        case{2,'gk'}
            model=2;
        case{3,'kmt'}
            model=3;
        case{'gal','cgmall'}
            type=varargin{i};
        otherwise
            if ~ischar(varargin{i})
                msg=num2str(varargin{i});
            else
                msg=varargin{i};
            end
            error('%s - unknown argument: %s',current_function().lower,msg);
    end
    i=i+1;
end


if strcmp(sim,'none')
    error('%s - missing sim.',current_function().lower);
elseif strcmp(type,'none')
    error('%s - missing type.',current_function().lower);
elseif model==0
    error('%s - missing model.',current_function().lower);
end


switch(lower(sim))
    case 'tng50'
        legendPanel=1;
    case 'tng100'
        legendPanel=2;
end

figPos=[ 800          42        1062         954];
colors=brewermap(8,'Set1');
cind=[2 5 3 1];
htag={'11-12','12-13','13-14','14-15'};
mtag={'8.3-9','9-10','10-11','11-12'};
modelTag={'BR', 'GK', 'KMT'};


%% in Stellar Mass Bins 

hf=figure('color','w','position',figPos);

t=tiledlayout(2,2);
h=[];

% switch(lower(type))
%     case 'gal'
%         yl=[-0.4 0.25];
%     case 'cgmall'
%         yl=[-2 3];
% end

xl=[0 3.6];
for k=1:4
    nexttile;
    hflag=true;
    h=[];
    for j=1:4
        
        profSt=profStr.byHostStar4(k).(type).hiProfN(model,j);
        msk=profSt.binCount>0;
        if any(msk) && sum(profSt.binCount)>10
            
            px=[profSt.xMedian(msk) fliplr(profSt.xMedian(msk))];
            py=[profSt.yQuarts(2,msk) fliplr(profSt.yQuarts(3,msk))];
            
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            if hflag;hold on; hflag=false;end
            h(end+1)=plot(profSt.xMedian(msk),profSt.yMedian(msk),'-',...
                'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        end
        
    end
    
    if k==legendPanel
        legend(h,'fontsize',14,'location','southeast','interpreter','latex');
    end
    
    
    xlim(xl);
    %ylim(yl);
    xl=xlim;
    yl=ylim;
    
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14,'box','on')
    
    
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('HI mass /stellar mass',16);
    
    
end

t.TileSpacing = 'compact';

sgtitle([sim ' ' type ' ' modelTag{model} ' - Median and 25-75 percentile'],...
    'Interpreter','latex','fontsize',18);

if printFlag
    fname=['hiMassN_radProf_' type '_med' '_stellarMass_' modelTag{model} '_' sim];
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end
%% Average

hf=figure('color','w','position',figPos);

h=[];
%yl=[-0.4 0.25];
% xl=[0 2.4];

t=tiledlayout(2,2);

for k=1:4
    nexttile;
    hflag=true;
    h=[];
    for j=1:4
        
        profSt=profStr.byHostStar4(k).(type).hiProfN(model,j);
        msk=profSt.binCount>0;
        if any(msk) && sum(profSt.binCount)>10
            
            px=[profSt.xMean(msk) fliplr(profSt.xMean(msk))];
            py=[profSt.yMean(msk)-0.5.*profSt.yStanDev(msk), ...
                fliplr(profSt.yMean(msk)+0.5.*profSt.yStanDev(msk))];
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            if hflag;hold on; hflag=false;end
            h(end+1)=plot(profSt.xMean(msk),profSt.yMean(msk),'-',...
                'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
            %          h(end+1)=errorbar(profSt.xMean,profSt.yMean,profSt.yStanDev./2,'-',...
            %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
            %
            
            %         h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassMedD(1,:,j),'-',...
            %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        end
        
    end
    
    if  k==legendPanel
        legend(h,'fontsize',14,'location','southeast','interpreter','latex');
    end
    
    
    xlim(xl);
    %ylim(yl);
    xl=xlim;
    yl=ylim;
    
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{*} =' mtag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14,'box','on')
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('HI mass /stellar mass',16);
    
    
end
t.TileSpacing = 'compact';

sgtitle([sim ' ' type ' ' modelTag{model} '- Mean and StdDev'],...
    'Interpreter','latex','fontsize',18);


if printFlag
    fname=['hiMassN_radProf_' type '_avg' '_stellarMass_' modelTag{model} '_' sim];
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end

%% in Host Mass Bins 
legendPanel=3;
hf=figure('color','w','position',figPos);

t=tiledlayout(2,2);
h=[];

switch(lower(type))
    case 'gal'
        yl=[-0.4 0.25];
    case 'cgmall'
        yl=[-2 3];
end

% xl=[0 2.4];
for k=1:4
    nexttile;
    hflag=true;
    h=[];
    for j=1:4
        
        profSt=profStr.byHostStar4(j).(type).hiProfN(model,k);
        msk=profSt.binCount>0;
        if any(msk) && sum(profSt.binCount)>10
            
            px=[profSt.xMedian(msk) fliplr(profSt.xMedian(msk))];
            py=[profSt.yQuarts(2,msk) fliplr(profSt.yQuarts(3,msk))];
            
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            if hflag;hold on; hflag=false;end
            h(end+1)=plot(profSt.xMedian(msk),profSt.yMedian(msk),'-',...
                'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
        end
        
    end
    
    if k==legendPanel
        legend(h,'fontsize',14,'location','southeast','interpreter','latex');
    end
    
    
    xlim(xl);
    %ylim(yl);
    xl=xlim;
    yl=ylim;
    
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14,'box','on')
    
    
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('HI mass /stellar mass',16);
    
    
end

t.TileSpacing = 'compact';

sgtitle([sim ' ' type ' ' modelTag{model} ' - Median and 25-75 percentile'],...
    'Interpreter','latex','fontsize',18);

if printFlag
    fname=['hiMassN_radProf_' type '_med' '_hostMass_' modelTag{model} '_' sim];
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end
%% Average

hf=figure('color','w','position',figPos);

h=[];
%yl=[-0.4 0.25];
% xl=[0 2.4];
% 
t=tiledlayout(2,2);

for k=1:4
    nexttile;
    hflag=true;
    h=[];
    for j=1:4
        
        profSt=profStr.byHostStar4(j).(type).hiProfN(model,k);
        msk=profSt.binCount>0;
        if any(msk) && sum(profSt.binCount)>10
            
            px=[profSt.xMean(msk) fliplr(profSt.xMean(msk))];
            py=[profSt.yMean(msk)-0.5.*profSt.yStanDev(msk), ...
                fliplr(profSt.yMean(msk)+0.5.*profSt.yStanDev(msk))];
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            if hflag;hold on; hflag=false;end
            h(end+1)=plot(profSt.xMean(msk),profSt.yMean(msk),'-',...
                'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',mtag{j} );
            %          h(end+1)=errorbar(profSt.xMean,profSt.yMean,profSt.yStanDev./2,'-',...
            %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
            %
            
            %         h(end+1)=plot(hprofs50.byHostStar4(k).xMed(:,j),hprofs50.byHostStar4(k).Gal.hiMassMedD(1,:,j),'-',...
            %             'linewidth',1.2,'color',colors(cind(j),:),'DisplayName',htag{j} );
        end
        
    end
    
    if  k==legendPanel
        legend(h,'fontsize',14,'location','southeast','interpreter','latex');
    end
    
    
    xlim(xl);
    %ylim(yl);
    xl=xlim;
    yl=ylim;
    
    grid
    text(xl(1)+0.05.*diff(xl),yl(1)+0.95.*diff(yl),...
        ['$ M_\mathrm{h} =' htag{k} '$'],'interpreter','latex',...
        'fontsize',16);
    
    set(gca,'fontsize',14,'box','on')
    xlabelmine('$r/R_\mathrm{vir}$',16);
    ylabelmine('HI mass /stellar mass',16);
    
    
end
t.TileSpacing = 'compact';

sgtitle([sim ' ' type ' ' modelTag{model} '- Mean and StdDev'],...
    'Interpreter','latex','fontsize',18);

if printFlag
    fname=['hiMassN_radProf_' type '_avg' '_hostMass_' modelTag{model} '_' sim];
    printout_fig(gcf,fname,'nopdf','v','printoutdir',outDir);
end
