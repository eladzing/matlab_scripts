%% extract profile from 2D Image. Use recursive division of cells into sub cells 

function prof=xray_proj_profile_fancy(boxx,imaj,binFac)

%boxx=1; %in Mpc
global NCELL
global hub 
%imaj=projI;

%% generate bins - should be larger than cell size
cellsize = (boxx/hub)/NCELL;
%binFac=10; 
binSize0=binFac*cellsize;

rmin=cellsize;
rmax=boxx/hub/2;

nBins=floor((rmax-rmin)/binSize0);
binSize=(rmax-rmin)/nBins;

binEdges=rmin:binSize:rmax;
bins=binEdges(1:end-1)+diff(binEdges)./2;

%% create  distance array 

cm=[0 0];
[meshX, meshY] = meshgrid(1:NCELL);

meshX = meshX - (NCELL+1)/2 -cm(1);
meshY = meshY - (NCELL+1)/2 -cm(2);
% Fix Units (to be in Mpc)
meshX = meshX * ((boxx/hub)/NCELL);
meshY = meshY * ((boxx/hub)/NCELL);


%rcube=sqrt(meshX.^2+meshY.^2) ; % r cube in Mpc

%rCells=reshape(rcube,[numel(rcube) 1]);
meshX=reshape(meshX,[numel(meshX) 1]);
meshY=reshape(meshY,[numel(meshY) 1]);
fluxCells=reshape(imaj,[numel(imaj) 1]);

%% generate vertice list 
vertsX=zeros(NCELL^2,4);
vertsY=zeros(NCELL^2,4);
vertsB=zeros(NCELL^2,4);
vertsX(:,1)=meshX-0.5*cellsize;
vertsX(:,2)=meshX+0.5*cellsize;
vertsX(:,3)=meshX+0.5*cellsize;
vertsX(:,4)=meshX-0.5*cellsize;

vertsY(:,1)=meshY-0.5*cellsize;
vertsY(:,2)=meshY-0.5*cellsize;
vertsY(:,3)=meshY+0.5*cellsize;
vertsY(:,4)=meshY+0.5*cellsize;

clear meshX meshY 
%% find the bin for each vert
eps=1e-10;
for i=1:4
    vertRad0=sqrt(vertsX(:,i).^2+vertsY(:,i).^2); 
    for k=1:nBins+1
        mask=abs(vertRad0-binEdges(k))<eps;
        vertsX(mask,i)=vertsX(mask,i)+sign(vertsX(mask,i)).*eps;
        vertsY(mask,i)=vertsY(mask,i)+sign(vertsY(mask,i)).*eps;
        
    end
        
    vertRad=sqrt(vertsX(:,i).^2+vertsY(:,i).^2);
    vB=ceil((vertRad-rmin)./(rmax-rmin).*nBins);
    %vB1=ceil((vertRad.*(1+eps)-rmin)./(rmax-rmin).*nBins);
    %vB=max(vB0,vB1);
    vB(vertRad==0)=1;
    vertsB(:,i)=vB;
end

binFlux=zeros(nBins,1);
binArea=zeros(nBins,1);
binCount=zeros(nBins,1);
areaTot=0;
area1=zeros(size(fluxCells));
area2=zeros(size(fluxCells));
wtbar=length(fluxCells);
hwb=waitbar(0);
stp=ceil(wtbar/100);
for j=1:length(fluxCells)
    if mod(j,stp)==0
        waitbar(j./wtbar,hwb);
    end
    
    xx=vertsX(j,:);
    yy=vertsY(j,:);
    bb=vertsB(j,:); 
    
    area1(j)=polygonArea(xx,yy);
    
    
%     if any(abs(xx+0.3571)<1e-4) && any(abs(yy)<1e-4) 
%         disp(bb);
%     end
%     
    if ~all(bb<=nBins & bb>=1)  %ignore cells which are out of bounds 
        continue 
    end
    
    % subdivide cells 
    if ~any(diff(bb)) % all vertices in the same bin 
        binFlux(bb(1))=binFlux(bb(1))+fluxCells(j).*area1(j);
        binArea(bb(1))=binArea(bb(1))+area1(j);
        binCount(bb(1))=binCount(bb(1))+1;
        areaTot=areaTot+area1(j);
        area2(j)=area1(j);
    else  % subdivide 
        celPnts=[];
        for i=1:length(xx)
            ip=i+1-length(xx).*(i+1>length(xx));
            pnts=[];
            pnts=subdivide_edge(pnts,[xx(i) xx(ip)],[yy(i) yy(ip)],[bb(i) bb(ip)],binEdges(2:end),0);
            celPnts=cat(1,celPnts,pnts);
        end
        
        for i=min(celPnts(:,3)):max(celPnts(:,3))
            pp=celPnts(celPnts(:,3)==i,1:2);
            area=polygonArea(pp(:,1),pp(:,2));
            binFlux(i)=binFlux(i)+fluxCells(j).*area;
            binArea(i)=binArea(i)+area;
            binCount(i)=binCount(i)+area./area1(j);
            areaTot=areaTot+area;
            area2(j)=area2(j)+area;
            
        end
    end
        
end
close(hwb)

prof.profile=binFlux./binArea;
prof.binEdges=binEdges;
prof.bins=bins;
prof.binArea=binArea;
prof.binCount=binCount;




end       
 


function points=subdivide_edge(points,xx,yy,bb,binRad,pInd)
%% recursive function to dividee the edge of a cell into segments
% according to the radial bins that edge crosses

% add initial point to the list
pInd=pInd+1;
points(pInd,1)=xx(1);
points(pInd,2)=yy(1);
points(pInd,3)=bb(1);

if diff(bb)~=0  % second point not in the same bin
    % determine direction of line w.r.t bin order
    upfac=diff(bb)/abs(diff(bb));
    
    if upfac>0
        rad=binRad(bb(1));
    else
        rad=binRad(bb(1)-1);
    end
    
    %find intersection point
    
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
    

function [xi,yi]=divide_edge(xx,yy,rad)
%% auxilary function for dividing a cell across several radial shells
% given 2 points (p1,p2) which define a cell edge ( orientation is
% important - line is from p1 to p2), find the intersection point between
% the line and the edge of a radial shell (it is already known that the
% points are in different bins

if (length(xx)~=2 || length(yy) ~=2)
    error('DIVIDE_EDGE: arguments must be of length 2');
end

%% find intersection points
if diff(xx)==0  %solve for vertical lines
    intrPnt(1).x=xx(1);
    intrPnt(2).x=xx(1);
    intrPnt(1).y=sqrt(rad^2-xx(1)^2);
    intrPnt(2).y=-sqrt(rad^2-xx(1)^2);
    
    
elseif diff(yy)==0 % solve for horizontal lines
    intrPnt(1).y=yy(1);
    intrPnt(2).y=yy(1);
    intrPnt(1).x=sqrt(rad^2-yy(1)^2);
    intrPnt(2).x=-sqrt(rad^2-yy(1)^2);
    
else
    
    %find line parameters
    
    a=diff(yy)/diff(xx);
    b=yy(1)*xx(2)-yy(2)*xx(1)/diff(xx);
    
    intrPnt=find_line_radius_intersection(a,b,rad);
end

%% choose the correct solution
if isstruct(intrPnt)
    if length(intrPnt)<=2
        found=false(1,2);
        for i=1:length(intrPnt)
            sgnX=(intrPnt(i).x-xx(1))/(xx(2)-intrPnt(i).x);
            sgnY=(intrPnt(i).y-yy(1))/(yy(2)-intrPnt(i).y);
            
            if sgnX >0 || sgnY>0
                found (i)=true;
            end
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
        error('DIVIDE_EDGE: neither intersection point found between endpoints')
    end
else
    error('DIVIDE_EDGE: both intersection points found between endpoints')
end


end

function point=find_line_radius_intersection(a,b,rad)
%% given a line defined by a & b (y=a*x+b), find the intersection
% with a circle of a given radius ( solve R^2=x^2+(a*x+b)^2

% avoid infinite descriminante
if ~isinf(a)
    descrim=a^2*b^2-(a^2+1)*(b^2-rad^2);
    fac=1/(a^2+1);
    
    eps=1e-8;
    
    if descrim>eps^2 %two real solutions
        deterSq=sqrt(descrim);
        point(1).x=fac*(-1*a*b+deterSq);
        point(2).x=fac*(-1*a*b-deterSq);
        
        point(1).y=a*point(1).x+b;
        point(2).y=a*point(2).x+b;
        
    elseif abs(descrim)<eps^2 %one real solution
        point(1).x=fac*(-1*a*b);
        point(1).y=a*point(1).x+b;
    elseif descrim<-eps^2
        point=Nan;
    end
else
    error('FIND_LINE_RADIUS_INTERSECTION: line is of infinite slope');
end
end
    

% %% divide image into bins 

% profile=zeros(nBins,1);
% for i=1:nBins 
%     mask=rcube>=binEdges(i) & rcube<binEdges(i+1);
%     
%     profile(i)=sum(imaj(mask))/(pi *( binEdges(i+1)^2-binEdges(i)^2))*cellsize^2;
% end
% 
