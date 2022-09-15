function hf = plot_demographics_2sims(maskJF,xVal,xBins,lineVal,lineBins,mask50,varargin)
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
xxGiven=[];

axFont=18;
legFont=20;
labFont=20;

legLoc={'NorthEast','northwest'};
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
            xxGiven=varargin{i};
        case{'leglab','legtag'}
            i=i+1;
            legTag=varargin{i};
        case{'legend'}
            i=i+1;
            fullLegend=varargin{i};
        case'axfont'
            i=i+1;
            axFont=varargin{i};
        case'legfont'
            i=i+1;
            legFont=varargin{i};
        case 'legloc'
            i=i+1;
            legLoc=varargin{i};
        case'labfont'
            i=i+1;
            labFont=varargin{i};
        otherwise
            error('%s - Unkown argument: %s',current_function().upper,varargin{i});
    end
    i=i+1;
end

%% find bin indices and x values 
xxBI=discretize(xVal,xBins);
lineBI=discretize(lineVal,lineBins);


%% prepare lines
flines=zeros(length(xBins)-1,length(lineBins)-1,2);
count=zeros(length(xBins)-1,length(lineBins)-1,2);
xx=zeros(length(xBins)-1,length(lineBins)-1,2);
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
            
            count(i,j,k)=sum(binMask);
            if sum(binMask)>0
                xx(i,j,k)=median(xVal(binMask));
                jf=maskJF(binMask);
                
                flines(i,j,k)=sum(jf)./length(jf);
                
                if sum(binMask)>1
                    bsci(i,j,:,k)=bootci(nboot,{fffun,jf});
                else
                    bsci(i,j,:,k)=jf;
                end
            else
                xx(i,j,k)=0.5*sum(xBins(i:i+1));
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

if logFlag
        xx(xx>0)=log10(xx(xx>0));
end


hf=myFigure();
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
        
        ind=find(count(:,j,k)>0,1,'first'):find(count(:,j,k)>0,1,'last');
        
        if isempty(xxGiven)
            xxx=xx(:,j,k)';
        else
            xxx=xxGiven;
        end
        
        
        if ~isempty(ind)
            
            px=[xxx(ind) fliplr(xxx(ind))];
            py=[bsci(ind,j,1,k)' fliplr(bsci(ind,j,2,k)')];
            patch(px,py,colors(cind(j),:),'facealpha',0.35,'edgecolor','none')
            if j==1
                hold on
                h2(k)=plot(xxx(ind)',flines(ind,j,k),ls(k),...
                    'color','k','linewidth',lw,...
                    'displayName',simTag(k));
                
            end
            h(end+1)=plot(xxx(ind)',flines(ind,j,k),ls(k),...
                'color',colors(cind(j),:),'linewidth',lw,...
                'displayName',tag);
            
            
        else
            if j==1
                hold on
                h2(k)=plot(xxx',flines(:,j,k),ls(k),...
                    'color','k','linewidth',lw,...
                    'displayName',simTag(k));
                
            end
            h(end+1)=plot(xxx',flines(:,j,k),ls(k),...
                'color',colors(cind(j),:),'linewidth',lw,...
                'displayName',tag);
        end
        
    end
end

nc=1;
% if length(h)/2<4
%     nc=1;
% end
    

hl1=legend(h(1:2:end),'interpreter','latex','fontsize',legFont,'location',legLoc{1},'NumColumns',nc);
grid
xlabelmine(xLab,labFont);
ylabelmine('JF Fraction',labFont);
set(gca,'fontsize',axFont,'box','on')

ah1=axes('position',get(gca,'position'),'visible','off');

hl2=legend(ah1,h2(1:2),'interpreter','latex','fontsize',legFont,'location',legLoc{2},'NumColumns',2);



%titlemine(simTag);





