function res = tree2d( xVal,yVal,blCorner,xSide,ySide,splitParam,minLevel )
%TREE2D - Create an adaptivly sized grid based based on the number of data
%points per bin 
% This function generates an adaptive sized grid based such that each
% recatgular bin contains no more than a given number of data points. 
% The process is recursive - each bin which contains too many points is
% split into 4 equally sized bins which conserve the aspect ratio of the 
% bin original (basically a tree-code) 
% Each split defines a refinement level - the orginal plane is level 0, the
% first split is level one and so on. A minimal refinement level can be set
% to avoid bins which are very large and sparse. 
% 
% Argumnets: 
%   xVal,yVal     - x-axis and y-axis values of data points
%   blCorner      - x,y coordinates of bottom-left corner of the entire grid
%   xSide, ySide  - size of data plane in the X & Y directions
%   splitParam    - Maximal no. of points to per bin before splitting the bin
%   minlevel      - The minimal refinement level  -
% 
% ouptut: A data structure containing the input arguments and: 
% treeBLC     - an array of the x,y coordinates for the bottom-left corner
%               of all grid cells
% treeLevel   - the refinement level of each grid cell 
% pointsBLC   - for each data point, the x,y coordinate (bottom-left corner)
%               of the grid cell it inhabits 
% pointsLevel - for each data point, the refinement level of the grid cell
%                it inhabits 


% helpful arrays 
ddx=[0 1 1 0];
ddy=[0 0 1 1];

if ~exist('minLevel','var')
    minLevel=0;
end

% initialize output 
treeBLC=[];  % record x,y values of bottom-left corners of bins
treeLevel=[]; % record the refinement level of bins 
outputBLC=zeros(2,length(xVal)); % rAssign each data point to a grid bin
outputLevel=zeros(size(xVal)); % Aecord the bin refinement level for each data point  


% begin recursion at level= 0 
level=0;

doCell(blCorner,level) % recursive function which does the splitting. 


% assifgn output to data structure 
res.blCorner=blCorner;
res.xSide=xSide;
res.ySide=ySide;
res.pointsBLC=outputBLC;
res.pointsLevel=outputLevel;
res.treeBLC=treeBLC;
res.treeLevel=treeLevel;



%% Auxillary function


    function doCell(blc,level)
    % Recursive function for dealing with a grid cell: 
    % First check how many points are in it then: 
    % - add it to the tree if the number of points are below the limit
    % - split it up (call the same function) if too many points 
    % grid cells with zero points are not recorded 
   
    % Check which points are found in the cell 
    inCellMask=inCell(blc,level);
    
    if sum(inCellMask)>0 
        if sum(inCellMask)<=splitParam && level>=minLevel
            
            % number of data points is below the limit & cell is above the 
            % minmal refinement level
            % This cell is done! 
                        
            outputBLC(1,inCellMask)=blc(1);
            outputBLC(2,inCellMask)=blc(2);
            outputLevel(inCellMask)=level;
            treeBLC(:,end+1)=blc;
            treeLevel(end+1)=level;
            
            
        else
            
            % split cell into 4 and repeat for each sub-cell
            
            blcNew(:,1)=blc(1)+ddx.*(xSide/2^(level+1));
            blcNew(:,2)=blc(2)+ddy.*(ySide/2^(level+1));
            
            for i=1:length(blcNew)
                doCell(blcNew(i,:),level+1)
            end
        end
    end
    
    
    end



    function mask= inCell(blc,level)
    % Auxilary function which returns a logical mask for which points are
    % found within a grid bin. 
    
    %% build polygon
    xPoly=blc(1)+ddx.*(xSide/2^level); %      [x0 x0+dx x0+dx x0];
    yPoly=blc(2)+ddy.*(ySide/2^level);%      [y0 y0 y0+dy y0+dy];
    
    
    
    mask=inpolygon(xVal,yVal,xPoly,yPoly);
    end

end
