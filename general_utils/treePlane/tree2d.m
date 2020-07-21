function res = tree2d( xVal,yVal,blCorner,xSide,ySide,splitParam,minLevel )
%UNTITLED Summary of this function goes here
%   xx,yy - x-axis and y-axis values of points
%   xSide, ySide  - size of data plane in the X & Y directions
%   splitParam - threshold no. of points to split


ddx=[0 1 1 0];
ddy=[0 0 1 1];

if ~exist('minLevel','var');
    minLevel=0;
end

outputBLC=zeros(2,length(xVal));
outputLevel=zeros(size(xVal));
treeBLC=[];
treeLevel=[];



level=0;

doCell(blCorner,level)


res.blCorner=blCorner;
res.xSide=xSide;
res.ySide=ySide;
res.pointsBLC=outputBLC;
res.pointsLevel=outputLevel;
res.treeBLC=treeBLC;
res.treeLevel=treeLevel;



%% Auxillary function


    function doCell(blc,level)
    
    
    inCellMask=inCell(blc,level);
    
    if sum(inCellMask)>0 
        if sum(inCellMask)<=splitParam && level>=minLevel
            
            outputBLC(1,inCellMask)=blc(1);
            outputBLC(2,inCellMask)=blc(2);
            outputLevel(inCellMask)=level;
            treeBLC(:,end+1)=blc;
            treeLevel(end+1)=level;
            
            
        else
            
            blcNew(:,1)=blc(1)+ddx.*(xSide/2^(level+1));
            blcNew(:,2)=blc(2)+ddy.*(ySide/2^(level+1));
            
            for i=1:length(blcNew)
                doCell(blcNew(i,:),level+1)
            end
        end
    end
    
    
    end




    function mask= inCell(blc,level)
    
    %% build polygon
    xPoly=blc(1)+ddx.*(xSide/2^level); %      [x0 x0+dx x0+dx x0];
    yPoly=blc(2)+ddy.*(ySide/2^level);%      [y0 y0 y0+dy y0+dy];
    
    
    
    mask=inpolygon(xVal,yVal,xPoly,yPoly);
    end

end
