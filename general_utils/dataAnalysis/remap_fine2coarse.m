
function [newMap, binsize, xxlim,yylim]= remap_fine2coarse(xx,yy,vv,cellSize0,varargin)
%% map a fine grid to a grid with larger cells
% This function maps a fine grid to a coarser one. Cells of the original
% grid which overlap more than one coarse cell divided into the appropriate
% contributions. The grids must both be cartesian, but not necassarily
% square.
% the result ('newMap') is a len(1) X len(2) array.


%The first index along the 3rd dimension is the
% 2D array for the sum of values which belong in each bin. Each value is
% weighted by the array wt which is one by defualt.
%
% The second index is an array containing the sum of weights for each bin.
% This allows one to find the average value in each bin.
%
% xx - 1st parameter array
% yy - 2nd parameter array
% vv - value to be binned
% varargin:
%  len - two vector of histogram bins for each dimension
%        defualt: 200 X 200
%  wt  - weight for each value. default: 1
%  xxlim,yylim - limiting values for the histogram. default is using
%                min/max of the arrays

len=[200 200]; % default dimensions of histogram
wt=ones(size(vv)); %default weight array

% set cell length
switch length(cellSize0)
    case 1
        cellSize=cellSize0.*[1 1];
    case 2
        cellSize=cellSize0;
    otherwise
        error('REMAP_FINE2COARSE - illegal cellsize')
end

xxlim=([min(xx) max(xx)]+0.5*[-1 1].*cellSize(1)).*1.000;
yylim=([min(yy) max(yy)]+0.5*[-1 1].*cellSize(2)).*1.000;

i=1;
while i<=length(varargin)
    switch varargin{i}
        case {'len','bins'}
            i=i+1;
            len=varargin{i};
            if length(len)==1
                len(2)=len;
            elseif (length(len)>2 || length(len)<1)
                error('histogram2d: size of len cna only be 1 or 2');
            end
        case {'xxlim','xlim'}
            i=i+1;
            xxlim=varargin{i};
        case {'yylim','ylim'}
            i=i+1;
            yylim=varargin{i};
        case{'wt','weight'}
            i=i+1;
            wt=varargin{i};
            if length(wt)==1
                wt=ones(size(xx)).*wt;
            end
        otherwise
            error('REMAP_FINE2COARSE - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

%% set the coarse grid
%av=[-1 1];
binsize=[diff(xxlim)/len(1) diff(yylim)/len(2)];
%xxlim=xxlim+0.5*binsize(1)*av;
%yylim=yylim+0.5*binsize(2)*av;
%len=len+1;
newMap=zeros(len(1),len(2),2);

% remove values which are out of bounds
ind=find((xx>=xxlim(1) & xx<=xxlim(2)) & (yy>=yylim(1) & yy<=yylim(2)));
xx=xx(ind);
yy=yy(ind);
vv=vv(ind);

vvw=vv.*wt;


% check to see that the original cellsize is smaller than the coarse grid
if any(cellSize>binsize)
    error('Cellsize must be smaller than binsize of coarse map')
end


%% assign fine grid cells to coarse grid cells according to cell centers.

indx=ceil((xx-xxlim(1))./binsize(1));
indy=ceil((yy-yylim(1))./binsize(2));

bottomLeft=[xxlim(1) yylim(1)];


%% run over cells to add contributions
hw=waitbar(0,'Remapping...');
barStep=floor(0.05.*length(xx));

for i=1:length(xx)
    
    ix=indx(i);
    iy=indy(i);
    
    % define vertices of the coarse cell
    v1=([ix-1 iy-1]).*binsize+bottomLeft;
    v2=v1+[1 0].*binsize;
    v3=v1+[1 1].*binsize;
    v4=v1+[0 1].*binsize;
    
    % find the distances of the fine cell center to the coarse cell sides
    d(1) = find_distance_to_line([xx(i) yy(i)],v1,v2);
    d(2) = find_distance_to_line([xx(i) yy(i)],v2,v3);
    d(3) = find_distance_to_line([xx(i) yy(i)],v3,v4);
    d(4) = find_distance_to_line([xx(i) yy(i)],v4,v1);
    cl(1:3)=cellSize(1);
    cl(2:4)=cellSize(2);
    % define area fractions
    
    fArea=0.5.*(1+min(2.*d./cl,ones(size(d))));
    
    % subdivide fine cell and assign contributions to coarse cells
    
    newMap(ix,iy,1)=newMap(ix,iy,1)+prod(fArea).*vvw(i);
    newMap(ix,iy,2)=newMap(ix,iy,2)+prod(fArea).*wt(i);
    
    
    if any(fArea<1.0)  % some sub-division required
        
        % ensure you don't spill over the edge of the coarse grid
        skipFlag=false(1,8);
        switch(ix)
            case 1
                skipFlag([1 4 6])=true;
            case len(1)
                skipFlag([3 5 8])=true;
        end
        switch(iy)
            case 1
                skipFlag([1 2 3])=true;
            case len(2)
                skipFlag([6 7 8])=true;
        end
        
        % distribute contributions
        for j=1:8
            if skipFlag(j)
                continue
            end
            
            % determine correct contribution and the cell addresses
            switch(j)
                case 1
                    fac=(1-fArea(1))*(1-fArea(4));
                    yInd=iy-1;
                    xInd=ix-1;
                    
                case 2
                    fac=(1-fArea(1))*fArea(2)*fArea(4);
                    yInd=iy-1;
                    xInd=ix;
                    
                case 3
                    fac=(1-fArea(1))*(1-fArea(2));
                    yInd=iy-1;
                    xInd=ix+1;
                    
                case 4
                    fac=(1-fArea(4))*fArea(1)*fArea(3);
                    yInd=iy;
                    xInd=ix-1;
                    
                case 5
                    fac=(1-fArea(2))*fArea(1)*fArea(3);
                    yInd=iy;
                    xInd=ix+1;
                    
                case 6
                    fac=(1-fArea(3))*(1-fArea(4));
                    yInd=iy+1;
                    xInd=ix-1;
                    
                case 7
                    fac=(1-fArea(3))*fArea(2)*fArea(4);
                    yInd=iy+1;
                    xInd=ix;
                    
                case 8
                    fac=(1-fArea(3))*(1-fArea(2));
                    yInd=iy+1;
                    xInd=ix+1;
            end
            %             if xInd<1 || yInd <1
            %                 pause
            %             end
            %
            %             if xInd>500 || yInd >500
            %                 pause
            %             end
            
            % add contribution to the correct cell
            newMap(xInd,yInd,1)=newMap(xInd,yInd,1)+fac.*vvw(i);
            newMap(xInd,yInd,2)=newMap(xInd,yInd,2)+fac.*wt(i);
        end
        
    end
    if mod(i,barStep)==0
        waitbar(i./length(xx),hw)
    end
end

close(hw)
if abs(sum(sum(newMap(:,:,2)))/sum(wt)-1)>1e-3
    error('remap_fine2coarse - problem with summing up - inconsistent value in weight')
end


end


