function res = findDistance(p1,p2,boxsize,dim)
%FINDDISTANCE find the distance between two (sets of) points, taking
%boundary conditions into account. If you want to find the distance between 
% a set of points and a center use centerObject
%   dim -  sets the dimensionality of the space - 3 by default. 

if ~exist('boxsize','var')
    boxsize=inf;
end

if ~exist('dim','var')
    dim=3;
end

if size(p1)~=size(p2)
    
    error('FINDDISTANCE - incompaticle sizes for p1,p2')
end

dm=find(size(p1)==dim);


dp=abs(p1-p2);

dp=min(dp,abs(dp-boxsize));

res=sqrt(sum(dp.^2,dm));


end

