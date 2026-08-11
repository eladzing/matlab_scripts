function res = find_k_nearest_neighbor_2D_vel(objectPositions,objectVelocities,kk,dVel,varargin)
%FIND_k_NEAREST_NEIBOR find the k-th nearest neighbor and the distance to it for a
% bunch of objects in the TNG box but in 2D projection + velocity space .
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

projDir=["x" "y" "z"];


queryPoints=objectPositions;
queryPointsVel=objectVelocities;

global LBox;
boxSize=LBox;
global cosmoStruct
global illUnits

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))

        case {'querypoints','points','qp'}
            i=i+1;
            queryPoints=varargin{i};
            i=i+1;
            queryPointsVel=varargin{i};
        case {'boxsize','box','boxlength'}
            i=i+1;
            boxSize=varargin{i};
        case {'projection','direction'}
            i=i+1;
            projDir=varargin{i};
        %case {}
        otherwise
            error('%s - Illegal argumnet: %s',current_function().upper,varargin{i})
    end
    i=i+1;
end

%% test to see that input arrays are legit 3D

sz1=size(objectPositions);
sz2=size(queryPoints);

if sz1(1)~=3
    error('%s - object positions must be 3D, of size 3 x N',...
        current_function().upper)
elseif sz2(1)~=3
    error('%s - query points must be 3D, of size 3 x M',...
        current_function().upper)
end

%% run over queryPoints and find k-th nearest neigbor

for i=1:length(queryPoints(1,:))

    %tic
    newPos=double(illustris.utils.centerObject(objectPositions,queryPoints(:,i),boxSize)); % in

    % calculate hubble flow

    vHub=cosmoStruct.hub.*100.*newPos.*illUnits.lengthUnit./1000; % in km/sec.

    totalVel=(objectVelocities+vHub)-queryPointsVel(:,i);

    for j=1:length(projDir)
        switch lower(projDir(j))

            case('x')
                projInd=1;
                inds=[2 3];
            case('y')
                projInd=2;
                inds=[1 3];
            case('z')
                projInd=3;
                inds=[1 2];
            otherwise
                error('%s - unknown projection direction: %s \n',current_function().upper,...
                    projDir(j));
        end
        vMask=find(abs(totalVel(projInd,:))<dVel);
        dist=sqrt(sum(newPos(inds,vMask).^2,1));
        [mn,ix]=mink(dist,kk+1);
        kkk=kk+double(mn(1)==0);
      
        res.distance(j,i)=mn(kkk);%.*illUnits.lengthUnit;
        %nearNeib.distanceNorm(i)=mn(2)./double(fofs.Group_R_Crit200(i));
        res.indx(j,i)=vMask(ix(kkk));
    end



    %nearNeib.m200c(i)=m200c(ix(2));
    %toc
end
end