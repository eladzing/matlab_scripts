function [ output_args ] = plot2dTreeMap(val,treeStruct,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


ddx=[0 1 1 0 0];
ddy=[0 0 1 1 0];
hf=[];

cmap=brewermap(256,'*Spectral');
valLim=[min(val) max(val);];
edgeFlag=false;

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case{'fig','figure','fighandle'}
            i=i+1;
            hf=varargin{i};
        case{'min','valmin'}
            i=i+1;
            valLim(1)=varargin{i};
        case{'max','valmax'}
            i=i+1;
            valLim(2)=varargin{i};
        case{'minmax','maxmin'}
            i=i+1;
            valLim=varargin{i};
        case{'colormap','cmap'}
            i=i+1;
            cmap=varargin{i};
        case{'edge'}
            edgeFlag=true;
            
        otherwise
            error('plot2dTreeMap - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

if isempty(hf)
    hf=figure;
else
    figure(hf);
end

%% scale colors

valScaled=(val-valLim(1))./diff(valLim);
valScaled(valScaled<0)=0;
%valScaled(isnan(valScaled))=0;
valScaled(valScaled>1)=1;
valScaled=ceil(valScaled.*(length(cmap)-1))+1;



% xPoly=treeStruct.blCorner(1)+ddx.*(treeStruct.xSide);
% yPoly=treeStruct.blCorner(2)+ddy.*(treeStruct.ySide);
%
% plot(xPoly,yPoly);




for i=1:length(treeStruct.treeLevel)
    
    
    if isnan(valScaled(i))
        continue
    end
    
    xPoly=treeStruct.treeBLC(1,i)+ddx.*(treeStruct.xSide/2^treeStruct.treeLevel(i));
    yPoly=treeStruct.treeBLC(2,i)+ddy.*(treeStruct.ySide/2^treeStruct.treeLevel(i));
    
    if ~edgeFlag
        patch(xPoly,yPoly,cmap(valScaled(i),:),'edgeColor','none');
    else
        patch(xPoly,yPoly,cmap(valScaled(i),:),'edgeColor','k','linewidth',0.1,'edgeAlpha',0.15);
        
    end
    
    
    if i==1
        hold on
    end
end

end

% res.blCorner=blCorner;
% res.xSide=xSide;
% res.ySide=ySide;
% res.pointsBLC=outputBLC;
% res.pointsLevel=outputLevel;
% res.treeBLC=treeBLC;
% res.treeLevel=treeLevel;
