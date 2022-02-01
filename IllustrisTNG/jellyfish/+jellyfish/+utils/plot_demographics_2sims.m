function hf = plot_demographics_2sims(maskJF,xxBI,xBins,lineBI,lineBins,mask50,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


nboot=1000;
cind=1:9;
lw=1.5;
colors=brewermap(9,'Set1');
%simMask=true(size(maskJF));
simTag=["TNG50" "TNG100"];
logFlag=false;
xLab='';
legTag='';
fullLegend='';
xx=[];

i=1;
while(i<=length(varargin))
    switch lower(varargin{i})
        
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
        case 'xx'
            i=i+1;
            xx=varargin{i};
        case{'leglab','legtag'}
            i=i+1;
            legTag=varargin{i};
        case{'legend'}
             i=i+1;
            fullLegend=varargin{i};
        otherwise
            error('%s - Unkown argument: %s',current_function().upper,varargin{i});
    end
    i=i+1;
end

%% prepare lines
flines=zeros(length(xBins)-1,length(lineBins)-1,2);
bsci=zeros(length(xBins)-1,length(lineBins)-1,2,2);
fffun = @(x)(sum(x)./length(x));

for k=1:2
    switch k
        case 1
            simMask=mask50;
        case 2
            simMask=~mask50;
    end
    for i=1:length(xBins)-1
        for j=1:length(lineBins)-1
                       
            binMask=xxBI==i & lineBI==j & simMask';
            
            if sum(binMask)>0
                jf=maskJF(binMask);
                
                flines(i,j,k)=sum(jf)./length(jf);
                
                if sum(binMask)>1
                    bsci(i,j,:,k)=bootci(nboot,{fffun,jf});
                else
                    bsci(i,j,:,k)=jf;
                end
            end
        end
    end
end



%% plot figures

if isempty(xx)
    xx=(xBins(1:end-1));
    if logFlag
        xx=log10(xx);
    end
end

px=[xx fliplr(xx)];


hf=figure('color','w');
h=[];
h2=[];

ls=["-","--"];
for j=1:length(lineBins)-1
    if isempty(fullLegend)
        tag=[legTag num2str(log10(lineBins(j))) '-' num2str(log10(lineBins(j+1)))  ];
    else
        tag=fullLegend(j);
    end
    
    for k=1:2
    
    py=[bsci(:,j,1,k)' fliplr(bsci(:,j,2,k)')];
    patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
    if j==1
        hold on
        h2(k)=plot(xx,flines(:,j,k),ls(k),...
        'color','k','linewidth',lw,...
        'displayName',simTag(k));
    end
    h(end+1)=plot(xx,flines(:,j,k),ls(k),...
        'color',colors(cind(j),:),'linewidth',lw,...
        'displayName',tag);
    
    
    end
end

hl1=legend(h(1:2:end),'interpreter','latex','fontsize',14,'location','northeast');
grid
xlabelmine(xLab);
ylabelmine('JF Fraction');
set(gca,'fontsize',14,'box','on')

ah1=axes('position',get(gca,'position'),'visible','off');

hl2=legend(ah1,h2(1:2),'interpreter','latex','fontsize',14,'location','northwest');



%titlemine(simTag);





