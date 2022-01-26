function hf = plot_demographics(maskJF,xxBI,xBins,lineBI,lineBins,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


nboot=500;
cind=1:9;
lw=1.5;
colors=brewermap(9,'Set1');
simMask=true(size(maskJF));
simTag='TNG100 \& TNG50';
logFlag=false;
xLab='';
legTag='';

i=1;
while(i<=length(varargin))
    switch lower(varargin{i})
        case 'tng100'
            simTag='TNG100';
            i=i+1;
            simMask=varargin{i};
        case 'log'
            logFlag=true;
        case 'nboot'
            i=i+1;
            nboot=varargin{i};
        case {'cind','palette','colors'}
            i=i+1;
            cind=varargin{i};
        case{'xlab','xlabel','label'}
            i=i+1;
            xLab=varargin{i};
        case{'legend','leglab','legtag'}
            i=i+1;
            legTag=varargin{i};
        otherwise
            error('%s - Unkown argument: %s',current_function.upper(),varargin{i});
    end
    i=i+1;
end

%% prepare lines
flines=zeros(length(xBins)-1,length(lineBins)-1);
bsci=zeros(length(xBins)-1,length(lineBins)-1,2);
fffun = @(x)(sum(x)./length(x));


for i=1:length(xBins)-1
    for j=1:length(lineBins)-1
        
        
        
        binMask=xxBI==i & lineBI==j & simMask';
        
        if sum(binMask)>0
            jf=maskJF(binMask);
            
            flines(i,j)=sum(jf)./length(jf);
            
            if sum(binMask)>1
                bsci(i,j,:)=bootci(nboot,{fffun,jf});
            else
                bsci(i,j,:)=jf;
            end
        end
    end
end



%% plot figures

xx=(xBins(1:end-1));
if logFlag
    xx=log10(xx);
end


px=[xx fliplr(xx)];


hf=figure('color','w');
h=[];

for j=1:length(lineBins)-1
    tag=[legTag num2str(log10(lineBins(j))) '-' num2str(log10(lineBins(j+1)))  ];
    
    py=[bsci(:,j,1)' fliplr(bsci(:,j,2)')];
    patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
    if j==1;hold on;end
    h(j)=plot(xx,flines(:,j),...
        'color',colors(cind(j),:),'linewidth',lw,...
        'displayName',tag);
end

legend(h,'interpreter','latex','fontsize',14,'location','best')
grid
xlabelmine(xLab);
ylabelmine('JF Fraction');
set(gca,'fontsize',14,'box','on')
titlemine(simTag);





