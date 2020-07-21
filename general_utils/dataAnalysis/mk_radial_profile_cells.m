function res = mk_radial_profile_cells(cellPos,cellRad,val,varargin)
%MK_RADIAL_PROFILE_CELLS Given a bunch of cells/particles, create a radial
% profile, averaged over shells
%   Detailed explanation goes here

cellPos=double(cellPos);
cellRad=double(cellRad);
val=double(val);


%% defaults

rmin=0;
rmax=inf;
nbins=100;
center=[0 0 0];
rRange=[];
wt=[];
profType='extensive';
logFlag=false;
inMask=true(size(cellRad));
%cellsize=1;

binFlag=false;  %% if true, use user supplied bins

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
        case{'nbins','nb'}
            i=i+1;
            nbins=varargin{i};
        case{'type'}
            i=i+1;
            profType=varargin{i};
        case{'wt','weight','weights'}
            i=i+1;
            wt=double(varargin{i});
            
            if size(wt)~=size(cellPos)
                error('MK_RADIAL_PROFILE_CELLS - weights must be same size as cell arrays');
            end
            
        case{'log'}
            logFlag=true;
        case{'mask'}
            i=i+1;
            inMask=varargin{i};
        otherwise
            error('mk_radial_profile_2d - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end

if isempty(profType)
    error('MK_RADIAL_PROFILE_CELLS - Profile type must be given extensive/instensive)');
end
%% set up stuff


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

%% prepare bins

if ~binFlag
    
    rp=rRange(1):diff(rRange)/nbins:rRange(2);
else
    rp=bins;
    nbins=length(rp)-1;
    %rRange=[rp(1) rp(end)];
end
rsh=0.5.*(rp(1:end-1)+rp(2:end));
binSize=diff(rp);
shellVol=(4*pi/3).*diff(rp.^3);
profile=zeros(size(rsh));
weights=zeros(size(profile));


%% remove out of bounds and assign bins
cLow=rCell-cellRad;
cHi=rCell+cellRad;

if logFlag
    cLow=log10(cLow);
    cHi=log10(cHi);
end

mask=(cHi>rp(1) & cLow<rp(end)) & inMask;

rCell=rCell(mask);
cellRad=cellRad(mask);
val=val(mask);
wt=wt(mask);
cLow=cLow(mask);
cHi =cHi(mask);

% find bin index for each cell, for the lower and higher edges, do not
% exceed profile edge
cellBins(1,:)=max(ceil((cLow-rp(1))./(rp(end)-rp(1)).*nbins),1);
cellBins(2,:)=min(ceil((cHi-rp(1))./(rp(end)-rp(1)).*nbins),nbins);

% thisdeals with cells whose lower bound exceeds below zero
if rp(1)==0 && ~logFlag
    negMask=rCell-cellRad<0;
    
else
    negMask=false(size(rCell));
end


%cellBins(1,cellBins(1,:)<1)=1;
%cellBins(2,cellBins(2,:)>nbins)=nbins;

%fracMask=diff(cellBins,1,1)~=0;

%% loop over cells, according to method  and add accordingly

switch profType
    case 'extensiveFancy'
        
        for i=1:length(rCell)
            ind=cellBins(1,i):cellBins(2,i);
            
            % get volujme fractions
            r1=rp(ind);
            r2=rp(ind+1);
            if logFlag
                r1=10.^r1;
                r2=10.^r2;
            end
            
            fracVol=partialSphereVolume(r1,r2,cellRad(i),rCell(i));
            
            if negMask(i)  % add contributions of cells extending below zero dump all in center.
                fracVol(1)=fracVol(1)+(1-sum(fracVol));
            end
            
            profile(ind)=profile(ind)+val(i).*fracVol;
            weights(ind)=weights(ind)+fracVol.*(4*pi/3*cellRad(i)^3);
        end
        
    case 'intensiveFancy'
        
        if isempty(wt)
            error('MK_RADIAL_PROFILE_CELLS - weights must be given for intensive parameters');
        end
        
        vv=val.*wt;
        
        
        for i=1:length(rCell)
            ind=cellBins(1,i):cellBins(2,i);
            
            % get volujme fractions r1=rp(ind);
            r1=rp(ind);
            r2=rp(ind+1);
            if logFlag
                r1=10.^r1;
                r2=10.^r2;
            end
            
            fracVol=partialSphereVolume(r1,r2,cellRad(i),rCell(i));
            
            if negMask(i)  % add contributions of cells extending below zero
                fracVol(1)=fracVol(1)+(1-sum(fracVol));
            end
            
            profile(ind)=profile(ind)+vv(i).*fracVol;
            weights(ind)=weights(ind)+wt(i).*fracVol;
        end
        
        profile=profile./weights;
        
    case 'average'
        
        
        for i=1:length(rCell)
            ind=cellBins(1,i):cellBins(2,i);
            
            profile(ind)=profile(ind)+val(i);
            weights(ind)=weights(ind)+1;
        end
        
        profile=profile./weights;
        
        
    case 'intensive'
        ind=min(max(ceil((rCell-rp(1))./(rp(end)-rp(1)).*nbins),1),nbins);
        vv=val.*wt;
        
        for i=1:length(rCell)
            %ind=cellBins(1,i):cellBins(2,i);
            
            profile(ind(i))=profile(ind(i))+vv(i);
            weights(ind(i))=weights(ind(i))+wt(i);
        end
        
        profile=profile./weights;
        
    case 'extensive'
        ind=min(max(ceil((rCell-rp(1))./(rp(end)-rp(1)).*nbins),1),nbins);
        
        
        for i=1:length(rCell)
            %ind=cellBins(1,i):cellBins(2,i);
            
            profile(ind(i))=profile(ind(i))+val(i);
            weights(ind(i))=weights(ind(i))+1;
        end
        
        
        
    otherwise
        error('MK_RADIAL_PROFILE_CELLS - Unknown profile type (should be extensive/instensive): %g',profType);
end

res.profile=profile;
res.rShell=rsh;
res.weights=weights;
res.shellVol=shellVol;
res.rBins=rp;
res.binSize=binSize;




% set r