
function [bird, binsize, xxlim,yylim]= histogram2d(xx,yy,vv,varargin)
%% 2D histogram
% This function creates a 2-D histogram - the result ('bird') is a
% len(1) X len(2) X 2 array. The first index along the 3rd dimension is the
% 2D array for the sum of values which belong in each bin. Each value is
% weighted by the array wt which is one by defualt.
%
% The second index is an array containing the sum of weights for each bin.
% This allows one to find the average value in each bin.
%
% xx - 1st parameter array
% yy - 2nd parameter array
% vv - value to be binned
% varargin:
%  len - two vector of histogram bins for each dimension
%        defualt: 200 X 200
%  wt  - weight for each value. default: 1
%  xxlim,yylim - limiting values for the histogram. default is using
%                min/max of the arrays

len=[200 200]; % default dimensions of histogram
wt=ones(size(vv)); %default weight array
xxlim=[min(xx) max(xx)];
yylim=[min(yy) max(yy)];

i=1;
while i<=length(varargin)
    switch varargin{i}
        case {'len','bins'}
            i=i+1;
            len=varargin{i};
            if length(len)==1
                len(2)=len;
            elseif (length(len)>2 || length(len)<1)
                error('histogram2d: size of len cna only be 1 or 2');
            end
        case {'xxlim','xlim'}
            i=i+1;
            xxlim=varargin{i};
        case {'yylim','ylim'}
            i=i+1;
            yylim=varargin{i};
        case{'wt','weight'}
            i=i+1;
            wt=varargin{i};
        otherwise
            error('histogram2d: Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

%% some preliminaries
if size(xx)~=size(yy)
    yy=xx';
end

av=[-1 1];
binsize=[diff(xxlim)/(len(2)-1) diff(yylim)/(len(1)-1)];

xxlimH=xxlim+0.5*binsize(1)*av;
yylimH=yylim+0.5*binsize(2)*av;
%len=len+1;

% remove values which are out of bounds
ind=find((xx>=xxlimH(1) & xx<=xxlimH(2)) & (yy>=yylimH(1) & yy<=yylimH(2)));
xx=xx(ind);
yy=yy(ind);
vv=vv(ind);
wt=wt(ind);

vvw=vv.*wt;
bird=zeros(len(1),len(2),2);

indx=ceil((xx-xxlimH(1))./binsize(1));
indy=ceil((yy-yylimH(1))./binsize(2));


%tic
for i=1:length(xx)   %% this is actually faster than going over the grid. 
    bird(indy(i),indx(i),1)=bird(indy(i),indx(i),1)+vvw(i);
    bird(indy(i),indx(i),2)=bird(indy(i),indx(i),2)+wt(i);
end
%toc

% % 
% % %% try new
%  tic
%  ind=cat(1,indy,indx);
%  bird2=accumarray(ind',vvw,len);
%  bird2wt=accumarray(ind',wt,len);
%  toc

end
