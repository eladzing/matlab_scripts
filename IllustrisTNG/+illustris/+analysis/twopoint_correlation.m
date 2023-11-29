function res = twopoint_correlation(pos,rprof,boxSize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('boxSize','var')
    global LBox
    boxSize=LBox;
end

ncount=zeros(1,length(rprof)-1);
for i=1:length(pos)

    newPos=illustris.utils.centerObject(pos,pos(:,i),boxSize);
    
    dist=sqrt(sum(newPos.^2,1));
    
    ncount=ncount+histcounts(dist,rprof); 
    %[ncount, rprof]=ncount+histcounts(dist); 
    
end

%% get shell volume 
shellVolume = illustris.analysis.shellVolume_2pcf(rprof,boxSize);

%% 
nObject=length(pos);

numberDens=nObject./boxSize^3;


res= ( ncount./shellVolume.shellVolume)./(nObject.*numberDens);
%res.rprof=rprof;


