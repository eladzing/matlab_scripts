function res=cell2grid(coords,vals,cellSize,varargin)
%CELL2GRID - map scattered cells into uiform grid
%            Map a set of values in cells (e.g. from an unstructured mesh
%            or particel data) to a uniform grid.size of grid is set by the
%            positions of the particles, assuming they are centered on the
%            origin. Optionally, if not centered boxSide is set by the minimal
%            and maximal positions of the cells/particles - (DANGER: this
%            does not take into account periodic boundry conditions!!)
%
%            coordinates, cell values and cell sizes (diameter not radius!) must be given.
%            if cellSize is a single number than it is expanded to the size of
%            the vals array.
%            varargin can set the grid size, values weights, and setting
%            whether the values are extensive or intensive, as well as
%            whether the coordinates are centered.
%
%            Once boxside is set by min/max of cell/particle positions, it
%            is increased by the equivilant of a set number of grid cells.
%            This ensures that large cells on the edge do not spill over.
%            the number of cell lengths added to the boxside is set by the
%            variable 'buffer' which is set to 10 by default but can be set
%            in varargin.


%% defaults and perliminatries
len=length(vals);

Ngrid=256;
wt=ones(size(vals));
cellNum=len;
mappingType='null';
centeredFlag=true;
boxSide=0;
buffer=10; % once boxside is set by  min/max of part

if length(cellSize)==1
    cellSize=cellSize.*ones(size(vals));
end


%% parse arguments
i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case{'side','boxside','box'}
            i=i+1;
            boxSide=varargin{i};
        case{'ngrid','ng','nc','ncell'}
            i=i+1;
            Ngrid=varargin{i};
            
            Ngrid=Ngrid +mod(Ngrid,2); % fix values to be only even
            
            
        case{'weights','weight','wt'}
            i=i+1;
            wt=varargin{i};
        case{'extensive','ext','sum'}
            mappingType='extensive';
        case{'intensive','int','avg','average'}
            mappingType='intensive';
        case{'max','maximal'}
            mappingType='maximal';
        case{'notcentered','nocentered'}
            centeredFlag=false;
            fprintf('CAUTION: DOES NOT ACCOUNT FOR PERIODIC BOUNDARIES!! \n');
        case{'buffer','bf','buff'}
            i=i+1;
            buffer=varargin{i};
        otherwise
            error('CELL2GRID - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

if strcmp(mappingType,'null')
    error('CELL2GRID - Must state mapping type: extensive, intensive or maximal');
end
if ~any(size(coords)==3)
    error('CELL2GRID - Coordinates are not 3D');
end

if length(size(coords))~=2
    error('CELL2GRID - Coordinates wrong shape');
end


%% find box side length

dim=find(size(coords)~=3);
if boxSide==0
    if centeredFlag   % particles are centered on origin.
        m1=abs(max(coords,[],dim));
        m2=abs(min(coords,[],dim));
        ll=max(cat(1,m1,m2));
        
        boxSide=ceil(2.*ll);
        %   boxSide=boxSide+buffer*boxSide/Ngrid; %increase boxsize so smoothed cells won't go out of bounds
        clear m1 m2 ll
        
    else
        error('Non-centered positions not yet possible')
        %         m1=max(max(coords,[],dim));
        %         m2=min(min(coords,[],dim));
        %
        %         boxSide=ceil(m1-m2);
        %         %    boxSide=boxSide+buffer*boxSide/Ngrid; %increase boxsize so smoothed cells won't go out of bounds
        %         clear m1 m2 ll
        
    end
end

boxSide=boxSide+buffer*boxSide/Ngrid; %increase boxsize so smoothed cells will not go out of bounds

%% find location of sim cells in grid
% the following assumes that the positions are centered on the origin -
% need to chnage for uncentered case.
indX=ceil((coords(1,:)./boxSide+0.5).*Ngrid);
indY=ceil((coords(2,:)./boxSide+0.5).*Ngrid);
indZ=ceil((coords(3,:)./boxSide+0.5).*Ngrid);


gcl=ceil(cellSize./(boxSide/Ngrid));     %"diameter" of sim cells in grid units
%if strcmp(mappingType,'extensive-gaussian')
%    gcl=3.*gcl;
%end

gcl=gcl+(1-mod(gcl,2)); % fix values to be only odd

gind=(gcl+1)/2;   % extract no. of cells below and above the center cell;


indXlo=indX-gind;indXhi=indX+gind;
indYlo=indY-gind;indYhi=indY+gind;
indZlo=indZ-gind;indZhi=indZ+gind;


% ignore cells whose center is out of bounds
skipFlag= indX > Ngrid | indX < 1 | indY>Ngrid | indY<1 | indZ>Ngrid | indZ<1 ;
%skipFlag= indXlo>Ngrid | indXhi < 1 | indYlo>Ngrid | indYhi<1 | indZlo>Ngrid | indZhi<1 ;

indXlo(indXlo<1)=1;
indYlo(indYlo<1)=1;
indZlo(indZlo<1)=1;

indXhi(indXhi>Ngrid)=Ngrid;
indYhi(indYhi>Ngrid)=Ngrid;
indZhi(indZhi>Ngrid)=Ngrid;


%% build grid
cube=zeros(Ngrid,Ngrid,Ngrid);

% define postions of grid cell vertices in real world
gl=boxSide/Ngrid;

xg=-0.5*boxSide:gl:0.5*boxSide;
yg=xg;zg=xg;
[vertY,vertX,vertZ]=meshgrid(xg,yg,zg);


%% set value and weights and fill up grid

switch(mappingType)
    case('extensive')
        
        for i=1:cellNum
            
            if skipFlag(i)
                continue;
            end
            
            % grid center indices
            xi=indXlo(i):indXhi(i);
            yi=indYlo(i):indYhi(i);
            zi=indZlo(i):indZhi(i);
            
            % grid vertices indices - +1 since there is always one more vertex
            xvi=indXlo(i):indXhi(i)+1;
            yvi=indYlo(i):indYhi(i)+1;
            zvi=indZlo(i):indZhi(i)+1;
            
            % find vertices enclosed within cell
            vertRad=sqrt((vertX(xvi,yvi,zvi)-coords(1,i)).^2 +...
                (vertY(xvi,yvi,zvi)-coords(2,i)).^2 +  ...
                (vertZ(xvi,yvi,zvi)-coords(3,i)).^2);
            
            
            vertMask=vertRad<=cellSize(i)/2;
            
            % for each grid cell, count how many of its 8 vertices is enclosed in cell
            vw=zeros(length(xi),length(yi),length(zi));
            
            for ii=0:1
                for jj=0:1
                    for kk=0:1
                        
                        vw=vw+vertMask(1+ii:end-1+ii,...
                            1+jj:end-1+jj,...
                            1+kk:end-1+kk);
                    end
                end
            end
            
            
            vw=vw./8; % fractional weight of each cell
            
            % fix cells which fit into a grid cell (central grid cell always =1)
            mid=floor(size(vw)./2)+1;
            
            vw(mid(1),mid(2),mid(3))=1;
            
            % find total ''volume''
            volTot=sum(sum(sum(vw)));
            rho=vals(i)./volTot;
            
            % add contribution.
            cube(xi,yi,zi)=...
                cube(xi,yi,zi)+rho.*vw;
            
            
            
        end
        
        
        wts=0;
    case('intensive')
        wts=zeros(Ngrid,Ngrid,Ngrid); % for intensive values the weights are divided
        
        %wt=wt./(gcl.^3);
        
        for i=1:cellNum
            
            if skipFlag(i)
                continue;
            end
            
            % grid center indices
            xi=indXlo(i):indXhi(i);
            yi=indYlo(i):indYhi(i);
            zi=indZlo(i):indZhi(i);
            
            % grid vertices indices - +1 since there is always one more vertex
            xvi=indXlo(i):indXhi(i)+1;
            yvi=indYlo(i):indYhi(i)+1;
            zvi=indZlo(i):indZhi(i)+1;
            
            % find vertices enclosed within cell
            vertRad=sqrt((vertX(xvi,yvi,zvi)-coords(1,i)).^2 +...
                (vertY(xvi,yvi,zvi)-coords(2,i)).^2 +  ...
                (vertZ(xvi,yvi,zvi)-coords(3,i)).^2);
            
            vertMask=vertRad<=cellSize(i)/2;
            % for each grid cell, count how many of its 8 vertices is enclosed in cell
            vw=zeros(length(xi),length(yi),length(zi));
            
            for ii=0:1
                for jj=0:1
                    for kk=0:1
                        
                        vw=vw+vertMask(1+ii:end-1+ii,...
                            1+jj:end-1+jj,...
                            1+kk:end-1+kk);
                    end
                end
            end
            
            vw=vw./8; % fractional weight of each cell
            
            % fix cells which fit into a grid cell (central grid cell always =1)
            mid=floor(size(vw)./2)+1;
            
            vw(mid(1),mid(2),mid(3))=1;
            
            % find total ''volume''
            volTot=sum(sum(sum(vw)));
            wtRho=wt(i).*vw/volTot;
            
            cube(xi,yi,zi)=...
                cube(xi,yi,zi)+vals(i).*wtRho;
            
            wts(xi,yi,zi)=...
                wts(xi,yi,zi)+wtRho;
            
            
        end
        cube=cube./wts;
        cube(wts==0)=0;
        
        
    case('extensive-gaussian')
        
        % define sub-grid around particle for gaussian calculation.
        
        % calculate contribution to each cell in sub-grid
        
        % add to master grid
        
        
        
    case('maximal')
        wts=0;
        
        for i=1:cellNum
            if skipFlag(i)
                continue;
            end
            
            xi=indXlo(i):indXhi(i);
            yi=indYlo(i):indYhi(i);
            zi=indZlo(i):indZhi(i);
            
            for ii=xi
                for jj=yi
                    for kk=zi
                        cube(ii,jj,kk)=max(cube(ii,jj,kk),vals(i));
                    end
                end
            end
        end
        
end


res.cube=cube;
res.Ngrid=Ngrid;
res.boxSide=boxSide;
res.cellVol=(boxSide/Ngrid)^3;
res.weights=wts;
res.mapping=mappingType;

end
