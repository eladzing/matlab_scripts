function [ output_args ] = plot2dTreeGrid(treeStruct,varargin)
%PLOT2DTREEGRID Plot an adaptively sized grid based on tree2d
%   Given an adaptively sized grid (genetrated by tree2d function), plot a
%   the grid cells on a plane.
%
% Arguments:
% treeStruct - adaptively sized grid cells data structure (generated with
%               tree2d function)
% additional, optional arguments set as pairs:
% 'fig' - figure handle

% set defaults
ddx=[0 1 1 0 0];
ddy=[0 0 1 1 0];
hf=[];

% parse arguments
i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case{'fig','figure','fighandle'}
            i=i+1;
            hf=varargin{i};
        otherwise
            error('plot2dTreeGrid - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

% open/assign figure
if isempty(hf)
    hf=figure;
else
    figure(hf);
end


% define level 0 cell
xPoly=treeStruct.blCorner(1)+ddx.*(treeStruct.xSide);
yPoly=treeStruct.blCorner(2)+ddy.*(treeStruct.ySide);

plot(xPoly,yPoly);

hold on

% plot individual cells
for i=1:length(treeStruct.treeLevel)
    
    xPoly=treeStruct.treeBLC(1,i)+ddx.*(treeStruct.xSide/2^treeStruct.treeLevel(i));
    yPoly=treeStruct.treeBLC(2,i)+ddy.*(treeStruct.ySide/2^treeStruct.treeLevel(i));
    
    plot(xPoly,yPoly);
    
end

end

% res.blCorner=blCorner;
% res.xSide=xSide;
% res.ySide=ySide;
% res.pointsBLC=outputBLC;
% res.pointsLevel=outputLevel;
% res.treeBLC=treeBLC;
% res.treeLevel=treeLevel;
