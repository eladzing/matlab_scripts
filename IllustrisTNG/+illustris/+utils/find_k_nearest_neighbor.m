function res = find_k_nearest_neighbor(objectPositions,kk,varargin)
%FIND_k_NEAREST_NEIBOR find the k-th nearest neighbor and the distance to it for a
%bunch of objects in the TNG box. 
%   Given a list of object positions find the k-th nearest neigbor to each 
%   query point. If query points aren't given it is assumed that the object
%   positions are also the query points. 
%   Because of the periodic boundary conditions a 'brute-force' method is
%   used (rather than a kd-tree approach), based on mink  - finction for
%   finding the first k minima. 
%   if the query point is part of the object list than the first minima 
%   will always be zero and thus the k+1 minima is what we are looking for.  
%   This function uses the centerObject function which is in simulation
%   units. If other units are to be used, the boxsize must also be given in
%   appropriate units. 
%   output is in the same units as entered. 

warning('%s - positions should be in simulation units. Otherwise, a boxsize must be given',current_function().upper);

%% set defaults and parse arguments 

queryPoints=objectPositions;

global LBox;
boxSize=LBox;


i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))

        case {'querypoints','points'}
            i=i+1;
            queryPoints=varargin{i};
        case {'boxsize','box','boxlength'}
            i=i+1;
            boxSize=varargin{i};
        otherwise
            error('%s - Illegal argumnet: %s',current_function().upper,varargin{i})
    end
    i=i+1;
end



%% run over queryPoints and find k-th nearest neigbor 

for i=1:length(queryPoints)
    
    %tic
    newPos=double(illustris.utils.centerObject(objectPositions,queryPoints(:,i),boxSize));
    dist=sqrt(sum(newPos.^2,1));
    [mn,ix]=mink(dist,kk+1);
    kkk=kk+double(mn(1)==0);
    res.distance(i)=mn(kkk);%.*illUnits.lengthUnit;
    %nearNeib.distanceNorm(i)=mn(2)./double(fofs.Group_R_Crit200(i));
    res.ID(i)=ix(kkk)-1;
    %nearNeib.m200c(i)=m200c(ix(2));
    %toc
end
end