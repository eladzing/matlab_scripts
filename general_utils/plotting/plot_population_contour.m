function res=plot_population_contour(xdata,ydata,varargin)
%PLOT_POPULATION_CONTOUR - plot a contour plot of a population 
%   Plot a population of objects in  the xdata-ydata plane, as a contour
%   plot, where each contour is labeled by the percent of the population it
%   encloses. An optional smoothing filter, created by the fspecial function can also be given.


%% set default
maskFlag=false;
filterFlag=false;
mask=true(size(xdata));
xxlim=[min(xdata) max(xdata)];
yylim=[min(ydata) max(ydata)];
lc=[0 10 20 30 40 50 60 70 80 90 100] ;%    0:10:100;%lc(end+1)=99.9;
brewerMap='*YlOrRd';
plotFlag=true;

%% parse arguments

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case 'mask'
            i=i+1;
            mask=varargin{i};
            
            maskFlag=true;
        case {'xlim','xxlim'}
            i=i+1;
            xxlim=varargin{i};
            
        case {'ylim','yylim'}
            i=i+1;
            yylim=varargin{i};
            
        case {'filt','filter','smooth'}
            i=i+1;
            filterFlag=true;
            filt=varargin{i}; % must be filter created with fspecial function.
            
        case {'color','colormap','brewermap'}
            i=i+1;
            brewerMap=varargin{i};
        case {'levels','levs'}
            i=i+1;
            lc=varargin{i};
        case {'noplot','noshow'}
            plotFlag=false;
        otherwise
            error('plot_population_contour - illegal argument: %s',varargin{i});
    end
    
    i=i+1;
end

if maskFlag  % need to  redo some stuff 
    
    xxlim=[min(xdata(mask)) max(xdata(mask))];
    yylim=[min(ydata(mask)) max(ydata(mask))];
    
end

%% create histogram, and restructure values. 
[data, binsize, xxlim,yylim]= histogram2d(xdata(mask),ydata(mask),ones(size(xdata(mask))),...
    'xxlim',xxlim,'yylim',yylim,'bins',256);
data=data(:,:,1);


% smooth data
if filterFlag
    data=imfilter(data,filt);
    res.filter=filt;
end    

data=data./sum(sum(data)).*100;

bb=zeros(size(data));

for i=1:size(bb,1)
    for j=1:size(bb,2)
        maskTmp=data>=data(i,j);
        bb(i,j)=sum(data(maskTmp));
    end
end
bb(bb==max(max(bb)))=100;




%% plot

%xx=xxlim(1)+0.5*binsize(1):binsize(1):xxlim(2)-0.5*binsize(1);
%yy=yylim(1)+0.5*binsize(2):binsize(2):yylim(2)-0.5*binsize(2);

xx=xxlim(1):binsize(1):xxlim(2);
yy=yylim(1):binsize(2):yylim(2);

if plotFlag
    
    figure
    map=brewermap(256,brewerMap);
    map(257,:)=[1 1 1];
    colormap(map);
    
    contour(xx,yy,bb,'ShowText','on','LineColor',[0 0 0],...
        'LevelList',lc,'Fill','on');
    set(gca,'Fontsize',12)
    
    grid
end

res.xx=xx;
res.yy=yy;
res.popContour=bb;


end
