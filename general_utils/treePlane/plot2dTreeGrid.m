function [ output_args ] = plot2dTreeGrid(treeStruct,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


ddx=[0 1 1 0 0];
ddy=[0 0 1 1 0];
hf=[];

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

if isempty(hf)
    hf=figure;
else
    figure(hf);
end


xPoly=treeStruct.blCorner(1)+ddx.*(treeStruct.xSide);
yPoly=treeStruct.blCorner(2)+ddy.*(treeStruct.ySide);

plot(xPoly,yPoly);

hold on

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
