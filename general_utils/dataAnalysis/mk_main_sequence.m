function res = mk_main_sequence(xParam,yParam,varargin)
%MK_MAIN_SEQUENCE generate a 'main-sequence' for 2 parameters
%   start with an initial relation and iteratively remove outliers (beyond
%   1 sigma)

%% set defaults
nBins=100;
xMin=min(xParam);
xMax=max(xParam);
binFlag=false;

fac=2;
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
        case{'fac','factor'}
            i=i+1;
            fac=varargin{i};
        otherwise
            error('%s - Illegal argument: %s',current_function().upper,varargin{i})
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

%% begin iterations 
itCnt=0;
itCntMax=50;
mask=true(size(yParam));
n0=sum(mask);
nn=0;

while(nn~=n0 )
    n0=sum(mask);
    ms=mk_meanMedian_bin(xParam(mask),yParam(mask),'bins',binEdges);
    
    for i=1:length(ms.xMean)
        inds=find(xParam>=binEdges(i) & xParam<binEdges(i+1));
        ylim=ms.yMean(i)+ms.yStanDev(i).*[-1 1].*fac;
        mm=yParam(inds)<ylim(1) | yParam(inds)>ylim(2);
        mask(inds(mm))=false;
    end
    nn=sum(mask);
    
    
    itCnt=itCnt+1;
    if itCnt>itCntMax
        error('%s - Does not converge',current_function().upper)
    end
    res.itr(itCnt).msX=ms.xMean;
    res.itr(itCnt).msY=ms.yMean;
    res.itr(itCnt).mask=mask;
end

res.msX=ms.xMean(ms.binCount>0);
res.msY=ms.yMean(ms.binCount>0);
res.mask=mask;