function res = mk_meanMedian_bin(xParam,yParam,varargin)
%MK_MEANMEDIAN_BIN find the mean and median in binned data
%   Given a 2-D data set defined by an x-paramter and y_parameter,
%   construct a mean/median relation in bins of the x-parameter

%% set defaults
nBins=100;
xMin=min(xParam);
xMax=max(xParam);
binFlag=false;



%% parse arguments
i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case{'nb','nbin','nbins'}
            i=i+1;
            nBins=varargin{i};
        case{'xmax'}
            i=i+1;
            xMax=varargin{i};
        case{'xmin'}
            i=i+1;
            xMin=varargin{i};
        case{'bins'}
            i=i+1;
            binEdges=varargin{i};
            binFlag=true;
            
        otherwise
            error('mk_meanMedian_bin - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end


%% setup bins
if ~binFlag
    db=(xMax-xMin)/nBins;
    binEdges=xMin:db:xMax;
end

db=diff(binEdges);
bins=binEdges(1:end-1)+0.5.*db;
nBins=length(bins);


%% find mean and median in bins
binYMean=zeros(1,nBins);
binYStd=zeros(1,nBins);
binYQuantile=zeros(5,nBins);
binXMean=zeros(1,nBins);
binXStd=zeros(1,nBins);
binXQuantile=zeros(5,nBins);
binCount=zeros(1,nBins);

binInd=discretize(xParam,binEdges);

for i=1:nBins
    mask=binInd==i;
    binYMean(i)=mean(yParam(mask));
    binYStd(i)=std(yParam(mask));
    binYQuantile(:,i)=quantile(yParam(mask),[0.1 0.25 0.5 0.75 0.9]);
    
    binXMean(i)=mean(xParam(mask));
    binXStd(i)=std(xParam(mask));
    binXQuantile(:,i)=quantile(xParam(mask),[0.1 0.25 0.5 0.75 0.9]);
    
    binCount(i)=sum(mask);
end

res.bins=bins;
res.binEdges=binEdges;
res.binCount=binCount;

res.yMean=binYMean;
res.yStanDev=binYStd;
res.yMedian=binYQuantile(3,:);
res.yQuarts=binYQuantile([1,2,4,5],:);

res.xMean=binXMean;
res.xStanDev=binXStd;
res.xMedian=binXQuantile(3,:);
res.xQuarts=binXQuantile([1,2,4,5],:);
res.quantiles=[0.1 0.25 0.75 0.9];

end

