function res = mk_radial_histogram2d_test(cellPos,cellRad,yy,vv,varargin)
%MK_RADIAL_HISTOGRAM2D Given a bunch of cells/particles, create a radial
% profile, averaged over shells
%   Detailed explanation goes here

cellPos=double(cellPos);
cellRad=double(cellRad);
yy=double(yy);
vv=double(vv);


%% defaults

len=[200 200]; % default dimensions of histogram
yylim=[min(yy) max(yy)];
wt=ones(size(vv)); %default weight array

rmin=0;
rmax=inf;

center=[0 0 0];
rRange=[];

%profType='extensive';
xLogFlag=false;
yLogFlag=false;
%cellsize=1;

binFlag=false;  %% if true, use user supplied bins
profFlag=false;
%% get arguments
i=1;
while i<=length(varargin)
    
    switch lower(varargin{i})
        case{'rmin'}
            i=i+1;
            rmin=varargin{i};
        case{'rmax'}
            i=i+1;
            rmax=varargin{i};
        case{'range','rrange','rads'}
            i=i+1;
            rRange=varargin{i};
        case {'center','cen'}
            i=i+1;
            center=varargin{i};
        case {'bins'}
            i=i+1;
            binFlag=true;
            bins=varargin{i};
        case{'len','nbins','nb'}
            i=i+1;
            len=varargin{i};
            if length(len)==1
                len(2)=len;
            elseif (length(len)>2 || length(len)<1)
                error('histogram2d: size of len cna only be 1 or 2');
            end
        case {'yylim','ylim'}
            i=i+1;
            yylim=varargin{i};
            
        case{'wt','weight','weights'}
            i=i+1;
            wt=varargin{i};
            
            if size(wt)~=size(cellPos)
                error('MK_RADIAL_PROFILE_CELLS - weights must be same size as cell arrays');
            end
            
        case{'xlog','logx'}
            xLogFlag=true;
        case{'ylog','logy'}    
            yLogFlag=true;
        case{'prof','profile'}
            % also calculate a profile which includes contributions from beyond yylim
            profFlag=true;
            
        otherwise
            error('mk_radial_profile_2d - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end


%% perliminaries

% arrange indices
sz=size(cellPos);
if sz(1)~=3
    cellPos=cellPos';
end

% generate cell radial positions
for k=1:3
    cellPos(k,:)=cellPos(k,:)-center(k);
end
rCell=sqrt(sum(cellPos.^2,1));


%% set up y axis limits and bins
av=[-1 1];
yBinSize=diff(yylim)/(len(1)-1);
yylimH=yylim+0.5*yBinSize*av;



%% set up radial limits and bins

% set up range
if isempty(rRange)
    if isinf(rmax)% set rmax if unset
        %set rmax to just above furthest cell
        rmax=1.01.*max(rCell+cellRad);
    end
    
    rRange=[rmin rmax];
else
    if length(rRange)~=2
        error('MK_RADIAL_PROFILE_CELLS - rRange must be a two value array');
    end
    
    if diff(rRange)<0
        tm=rRange(1);
        rRange(1)=rRange(2);
        rRange(2)=tm;
        clear tm
    end
end

%% prepare radial  bins

if ~binFlag
    rp=rRange(1):diff(rRange)/len(2):rRange(2);
else
    rp=bins;
    rRange=[rp(1) rp(end)];
end
rsh=0.5.*(rp(1:end-1)+rp(2:end));
xBinSize=diff(rp);
shellVol=(4*pi/3).*diff(rp.^3);

if profFlag
    profile=zeros(size(rsh));
    profWeights=zeros(size(profile));
end

nbins=len(2);

if length(rp)~=len(2)+1
    error('MK_RADIAL_HISTOGRAM2D - rp and len(2) are not the same');
end
%% remove out of bounds and assign bins
cLow=rCell-cellRad;
cHi=rCell+cellRad;

if xLogFlag
    cLow=log10(cLow);
    cHi=log10(cHi);
end


%% choose loop according to profFlag

if ~profFlag
    
    mask=(cHi>rp(1) & cLow<rp(end)) & (yy>=yylimH(1) & yy<=yylimH(2));
    
    
    
    rCell=rCell(mask);
    cellRad=cellRad(mask);
    yy=yy(mask);
    cLow=cLow(mask);
    cHi =cHi(mask);
    wt=wt(mask);
    vv=vv(mask);
    
    % find bin index for each cell, for the lower and higher edges, do not
    % exceed profile edge
    cellBins(1,:)=max(ceil((cLow-rp(1))./(rp(end)-rp(1)).*nbins),1);
    cellBins(2,:)=min(ceil((cHi-rp(1))./(rp(end)-rp(1)).*nbins),nbins);
    
    % this deals with cells whose lower bound exceeds below zero
    if rp(1)==0 && ~xLogFlag
        negMask=rCell-cellRad<0;
        
    else
        negMask=false(size(rCell));
    end
    
    
    %% loop over cells, according to method  and add accordingly
    indy=ceil((yy-yylimH(1))./yBinSize);
    
    bird=zeros(len(1),len(2));
    weights=zeros(len(1),len(2));
    
    vvw=vv.*wt;
    
    
    for i=1:length(rCell)
        xInd=cellBins(1,i):cellBins(2,i);
        
        % get volujme fractions r1=rp(ind);
        try
        r1=rp(xInd);
        r2=rp(xInd+1);
        if xLogFlag
            r1=10.^r1;
            r2=10.^r2;
        end
        
        catch
            pause
        end
            
        
        fracVol=partialSphereVolume(r1,r2,cellRad(i),rCell(i));
        
        if negMask(i)  % add contributions of cells extending below zero
            fracVol(1)=fracVol(1)+(1-sum(fracVol));
        end
        
        
        bird(indy(i),xInd)= bird(indy(i),xInd)+vvw(i).*fracVol;
        weights(indy(i),xInd)= weights(indy(i),xInd)+wt(i).*fracVol;
        
    end

        
else
    profMask=cHi>rp(1) & cLow<rp(end);
        
    rCell=rCell(profMask);
    cellRad=cellRad(profMask);
    yy=yy(profMask);
    cLow=cLow(profMask);
    cHi =cHi(profMask);
    wt=wt(profMask);
    vv=vv(profMask);
    birdMask=yy>=yylimH(1) & yy<=yylimH(2);
    
    % find bin index for each cell, for the lower and higher edges, do not
    % exceed profile edge
    cellBins(1,:)=max(ceil((cLow-rp(1))./(rp(end)-rp(1)).*nbins),1);
    cellBins(2,:)=min(ceil((cHi-rp(1))./(rp(end)-rp(1)).*nbins),nbins);
    
    % this deals with cells whose lower bound exceeds below zero
    if rp(1)==0 && ~xLogFlag
        negMask=rCell-cellRad<0;
        
    else
        negMask=false(size(rCell));
    end
    
    
    %% loop over cells, according to method  and add accordingly
    indy=ceil((yy-yylimH(1))./yBinSize);
    
    bird=zeros(len(1),len(2));
    weights=zeros(len(1),len(2));
    
    vvw=vv.*wt;
    if yLogFlag
        yy=10.^(yy);
    end
    yyv=yy.*vv;
    
    
    for i=1:length(rCell)
        xInd=cellBins(1,i):cellBins(2,i);
        
        % get volujme fractions r1=rp(ind);
        r1=rp(xInd);
        r2=rp(xInd+1);
        if xLogFlag
            r1=10.^r1;
            r2=10.^r2;
        end
        
        fracVol=partialSphereVolume(r1,r2,cellRad(i),rCell(i));
        
        if negMask(i)  % add contributions of cells extending below zero
            fracVol(1)=fracVol(1)+(1-sum(fracVol));
        end
        
        profile(xInd)=profile(xInd)+yyv(i).*fracVol;
        profWeights(xInd)=profWeights(xInd)+vv(i).*fracVol;
        
        if birdMask(i)
            try
                bird(indy(i),xInd)= bird(indy(i),xInd)+vvw(i).*fracVol;
            catch
                pause
            end
            weights(indy(i),xInd)= weights(indy(i),xInd)+wt(i).*fracVol;
        end
        
    end
    profile=profile./profWeights;
        
    
    
    
    
    
    
    
end

%% fill result structure
res.bird=bird;
res.rShell=rsh;
res.weights=weights;
res.rBins=rp;
res.xBinSize=xBinSize;
res.yBinSize=yBinSize;
res.rRange=rRange;
res.yylim=yylim;

if profFlag

    res.profile=profile;
    res.profWeights=profWeights;
    res.shellVol=shellVol;
end
    


