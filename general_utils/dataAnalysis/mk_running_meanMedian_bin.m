function res = mk_running_meanMedian_bin(xParam,yParam,varargin)
%MK_RUNNING_MEANMEDIAN_BIN find the mean and median in binned data
%   Given a 2-D data set defined by an x-paramter and y_parameter,
%   construct a mean/median relation in bins of the x-parameter

%% set defaults
nBins=1000;
xMin=min(xParam);
xMax=max(xParam);
%binFlag=false;

binSize=[];


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
        case{'dx','binsize','db'}
            i=i+1;
            binSize=varargin{i};
%         case{'bins'}
%             i=i+1;
%             binEdges=varargin{i};
%             binFlag=true;
%  
        otherwise
            error('mk_meanMedian_bin - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end


%% setup bins
if isempty(binSize)
    binSize=25*(xMax-xMin)/nBins;
end


bins=xMin:(xMax-xMin)/(nBins-1):xMax;

binEdges(1,:)=bins-0.5*binSize;
binEdges(2,:)=bins+0.5*binSize;

% db=diff(binEdges);
% bins=binEdges(1:end-1)+0.5.*db;
% nBins=length(bins);


%% find mean and median in bins
binYMean=zeros(1,nBins);
binYStd=zeros(1,nBins);
binYQuantile=zeros(3,nBins);
binXMean=zeros(1,nBins);
binXStd=zeros(1,nBins);
binXQuantile=zeros(3,nBins);

for i=1:nBins
    mask=xParam>=binEdges(1,i) & xParam<binEdges(2,i);
    binYMean(i)=mean(yParam(mask));
    binYStd(i)=std(yParam(mask));
    binYQuantile(1:3,i)=quantile(yParam(mask),[0.25 0.5 0.75]);
    
    binXMean(i)=mean(xParam(mask));
    binXStd(i)=std(xParam(mask));
    binXQuantile(1:3,i)=quantile(xParam(mask),[0.25 0.5 0.75]);
end

res.binSize=binSize;
res.bins=bins;
res.binEdges=binEdges;
res.yMean=binYMean;
res.yStanDev=binYStd;
res.yMedian=binYQuantile(2,:);
res.yQuarts=binYQuantile([1,3],:);

res.xMean=binXMean;
res.xStanDev=binXStd;
res.xMedian=binXQuantile(2,:);
res.xQuarts=binXQuantile([1,3],:);


end

