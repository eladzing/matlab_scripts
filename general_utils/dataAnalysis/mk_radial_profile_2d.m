
function res = mk_radial_profile_2d(map,varargin)
%MK_RADIAL_PROFILE_2D Given a 2D map of something, create a radial profile, averaged over shells
%   Detailed explanation goes here
%   mapSide is a 3 value vector the bottom left corner and mapsid


%% defaults




rmin=0;
rmax=inf;
nbins=100;
center=[0 0];
rRange=[];
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
        case{'bottomleft'}
            i=i+1;
            bottomLeft=varargin{i};
        case{'mapsize','mapside','side','sides'}
            i=i+1;
            mapSide=varargin{i};
            if length(mapSide)==1
                mapSide=mapSide.*[1 1];
            end
        otherwise
            error('mk_radial_profile_2d - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end

if (isempty(bottomLeft) || isempty(mapSide))
    error('mk_radial_profile_2d - Must supply map geometry: bottomLeft and mapSide');
end

%% generate position grids for vertices

mapSize=fliplr(size(map));

cellsize=mapSide./mapSize;

xv=(0:mapSize(1)).*cellsize(1)+bottomLeft(1);
yv=(0:mapSize(2)).*cellsize(2)+bottomLeft(2);

[xVerts, yVerts]=meshgrid(xv,yv);

% generate rcube
rVerts=sqrt(xVerts.^2+yVerts.^2);

if isempty(rRange)
    if isinf(rmax)% set rmax if unset
        %set rmax to maximal radius which is fully represented in the map
        mids=(size(rVerts)+1)./2;
        
        rmax=min(rVerts(mids(1),1),rVerts(1,mids(2)));
        
    end
    
    rRange=[rmin rmax];
else
    if length(rRange)~=2
        error('MK_RADIAL_PROFILE_2D - rRange must be a two value array');
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
end
rsh=0.5.*(rp(1:end-1)+rp(2:end));
binSize=diff(rp);

%% assign a bin to each vertex of the cells
%vertBins=zeros(size(rVerts));
vertBins=ceil((rVerts-rRange(1))./diff(rRange).*nbins);

%for k=1:length(rp)-1
%    mask=rVerts>=rp(k) & rVerts<rp(k+1);
%    vertBins(mask)=k;
%end

%% analyze
binMass=zeros(size(rsh));
binVol=zeros(size(rsh));

wtbar=mapSize(1);
hwb=waitbar(0);
stp=ceil(wtbar/100);

% run over cells
for i=1:mapSize(1)
    if mod(i,stp)
        waitbar(i./wtbar,hwb,sprintf('analyzing cells...'));
    end
    % get cell vertices
    ivx=[i i+1 i+1 i];
    xx=xVerts(1,ivx);
    %     fprintf('i=%i \n',i)
    for j=1:mapSize(2)
        
        
        
        % get cell vertices
        ivy=[j j j+1 j+1];
        yy=yVerts(ivy,1)';
        
        vbins=zeros(size(xx));
        for ib=1:length(ivx)
            
            vbins(ib)=vertBins(ivy(ib),ivx(ib));
        end
        % check to see if cell is out of bounds
        if all(vbins>nbins) || all(vbins)<=0
            continue
        end
        
        % check to see if entire cell is in a bin
        cellVol=polygonArea(xx,yy);
        if ~any(diff(vbins))
            binMass(vbins(1))=binMass(vbins(1))+map(j,i);
            binVol(vbins(1))=binVol(vbins(1))+cellVol; % calculate cell volume;
            
        else %% subdivide cells
            celPnts=[];
            
            for ii=1:length(xx)
                ip=ii+1-length(xx).*(ii+1>length(xx));
                Pnts=[];
                Pnts=subdivide_edge(Pnts,[xx(ii) xx(ip)],[yy(ii) yy(ip)],[vbins(ii) vbins(ip)],rp(2:end),0);
                celPnts=cat(1,celPnts,Pnts);
            end
            
            %extract volumes
            for ii=min(celPnts(:,3)):max(celPnts(:,3))
                if ii<1 || ii>nbins
                    continue;
                end
                pp=celPnts(celPnts(:,3)==ii,1:2);
                vol=polygonArea(pp(:,1),pp(:,2)); % calculate cell volume
                
                binMass(ii)=binMass(ii)+map(j,i).*(vol/cellVol);
                binVol(ii)=binVol(ii)+vol;
            end
        end
        
    end
end
close(hwb);


res.rp=rp;
res.rsh=rsh;
res.prof=binMass;
res.shellVol=binVol;

end


%% auxilary functions

function points=subdivide_edge(points,xx,yy,bb,binRad,pInd)
%% recursive function to divide the edge of a cell into segments according
% to the3 radial bins that edge crosses

% add initial point to the list
pInd=pInd+1;
points(pInd,1)=xx(1);
points(pInd,2)=yy(1);
points(pInd,3)=bb(1);



% see if second point not in the same bin OR if both points are beyond the
% maximal radius
if ~( diff(bb)==0 || all(bb>length(binRad)))
    % determine direction of line w.r.t bin order
    upfac=diff(bb)/abs(diff(bb));
    if upfac>0
        ind=bb(1);
        %rad=binRad(bb(1));
        
    else
        ind=bb(1)-1;
        ind=min(ind,length(binRad));
        %         rad=binRad(bb(1)-1);
        
    end
    rad=binRad(ind);
    % find intersection point
    [xn,yn]=divide_edge(xx,yy,rad);
    pInd=pInd+1;
    points(pInd,1)=xn;
    points(pInd,2)=yn;
    points(pInd,3)=bb(1);
    
    points=subdivide_edge(points,[xn xx(2)],[yn yy(2)],[bb(1)+upfac bb(2)],binRad,pInd);
    
else
    
    pInd=pInd+1;
    points(pInd,1)=xx(2);
    points(pInd,2)=yy(2);
    points(pInd,3)=bb(2);
    
end
end




function  [xi,yi]=divide_edge(xx,yy,rad)
%% auxilary function for dividing a cell across several shells
% given 2 vertices of a cell edge, p1 & p2, such that the line goes from p1
% to p2 (orientation is important), find the intersection b/w the line and
% the edge of a radial shell given by rad. It should already be established
% that the two points belong in different bins.

if(length(xx)~=2 || length(yy)~=2)
    error('DIVIDE_EDGE: arguments must be of length 2');
end



%% find intersection points
eps=1e-10;
test=abs(sqrt(xx.^2+yy.^2)./rad-1)<eps;
if any(test)  % see if endpoints intersect radius
    xi=xx(test);
    yi=yy(test);
    
    if sum(test)>1  % both points are the same
        if diff(xi)<eps && diff(yi)< eps
            xi=xi(1);
            yi=yi(1);
        else
            error('DIVIDE_EDGE: Strange error that probably will not happen');
        end
    end
    
else
    
    % find intersection between endpoints
    
    if ~diff(xx)==0 %avoid straight vertical lines
        
        % find line parameters
        [a,b]=make_line_function([xx(1) yy(1)],[xx(2) yy(2)]);
        
        % find intersection poins
        intrPnt=find_line_radius_intersection(a,b,rad);
    else
        
        intrPnt(1).x=xx(1);
        intrPnt(2).x= intrPnt(1).x;
        intrPnt(1).y=sqrt(rad^2-xx(1)^2);
        intrPnt(2).y=-1.*intrPnt(1).y;
        
    end
    %% choose the correct solutions
    if isstruct(intrPnt)
        
        if length(intrPnt)<=2
            found=false(1,2);
            for i=1:length(intrPnt)
                
                dotProduct=(intrPnt(i).x-xx(1))*(intrPnt(i).x-xx(2))+(intrPnt(i).y-yy(1))*(intrPnt(i).y-yy(2));
                
                if dotProduct<0
                    found(i)=true;
                elseif dotProduct==0
                    error('Was not expecting this...')
                end
                
                %                 sgnX=(intrPnt(i).x-xx(1))/(xx(2)-intrPnt(i).x);
                %                 sgnY=(intrPnt(i).y-yy(1))/(yy(2)-intrPnt(i).y);
                %                 if isinf(sgnY) && sgnY>0
                %                     sgnY=(yy(2)-intrPnt(i).y)/(intrPnt(i).y-yy(1));
                %                 end
                %                 if sgnX>=0 || sgnY>=0
                %
                %                 end
            end
        elseif length(intrPnt)>2
            error('DIVIDE_EDGE: more than 2 intersection candidates found');
        end
    else
        error('DIVIDE_EDGE: no intersection point found');
    end
    
    if ~all(found)
        if any(found)
            intrPnt=intrPnt(found);
            xi=intrPnt.x;
            yi=intrPnt.y;
        else
            error('DIVIDE_EDGE: neither intersection point found between  endpoints');
        end
    else
        if intrPnt(1).y==intrPnt(2).y &&  intrPnt(1).x==intrPnt(2).x
            xi=intrPnt(1).x;
            yi=intrPnt(1).y;
        else
            
            
            error('DIVIDE_EDGE: both intersection points found to be between the initial endpoints');
        end
    end
    
end
end


function point=find_line_radius_intersection(a,b,rad,cen)
%% FIND_LINE_RADIUS_INTERSECTION Given a line defined by parameters a,b (y=a*x+b), find the intersection with a circle of a given radius.
% This entails solving the equation R^2=(x-xc)^2+(a*x+b-yc)^2

if ~exist('cen','var')
    cen=[0 0];
end

qa=a^2+1;
qb=2*(a*(b-cen(2))-cen(1));
qc=(b-cen(2))^2+cen(1)^2-rad^2;

%avoid infinite descriminant
if ~isinf(a)
    
    descrim=qb^2-4*qa*qc;
    
    eps=1e-9;
    
    if descrim>eps^2  % two real solutions
        point(1).x=(-1*qb+sqrt(descrim))/(2*qa);
        point(2).x=(-1*qb-sqrt(descrim))/(2*qa);
        
        point(1).y=a.*point(1).x+b;
        point(2).y=a.*point(2).x+b;
        
    elseif abs(descrim)<eps^2 % one real solution
        point(1).x=(-1*qb)/(2*qa);
        point(1).y=a.*point(1).x+b;
    elseif descrim<-eps^2 % no real solutions
        point=Nan;
    end
else
    error( 'FIND_LINE_RADIUS_INTERSECTION: line is of infinite slope')
end
end


