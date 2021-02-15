
function [bird, binsize, xxlimH]= histogram1d(xx,vv,varargin)
%% 2D histogram
% This function creates a 2-D histogram - the result ('bird') is a
% len(1) X 2 array. The first index along the 2nd dimension is the
% 1D array for the sum of values which belong in each bin. Each value is
% weighted by the array wt which is one by defualt.
%
% The second index is an array containing the sum of weights for each bin.
% This allows one to find the average value in each bin.
%
% xx - parameter array
% vv - value to be binned
% varargin:
%  len -  histogram bins for defualt: 200
%  wt  - weight for each value. default: 1
%  xxlim - limiting values for the histogram. default is using
%                min/max of the arrays

len=200; % default dimensions of histogram
wt=ones(size(vv)); %default weight array
xxlim=[min(xx) max(xx)];

if diff(xxlim)==0
    xxlim=xxlim.*[0.99 1.01];
end

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case {'len','bins'}
            i=i+1;
            len=varargin{i};
            
        case {'xxlim','xlim'}
            i=i+1;
            xxlim=varargin{i};
        
        case{'wt','weight'}
            i=i+1;
            wt=varargin{i};
        otherwise
            error('histogram1d: Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

%% some preliminaries

av=[-1 1];
binsize=diff(xxlim)/(len-1);

xxlimH=xxlim+0.5*binsize(1)*av;

% remove values which are out of bounds
ind=find(xx>=xxlimH(1) & xx<=xxlimH(2));

xx=xx(ind);
vv=vv(ind);
wt=wt(ind);

vvw=vv.*wt;
bird=zeros(len,2);

indx=discretize(xx,xxlimH(1):binsize:xxlimH(2)); ceil((xx-xxlimH(1))./binsize);


%tic
for i=1:length(xx)   %% this is actually faster than going over the grid. 
    bird(indx(i),1)=bird(indx(i),1)+vvw(i);
    bird(indx(i),2)=bird(indx(i),2)+wt(i);
end
%toc

 
%% try new - Actually loop is faster than accumarray
% tic
% 
% bird2=accumarray(indx',vvw);
% bird2wt=accumarray(indx',wt);
% toc

% pause
end
