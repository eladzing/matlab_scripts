function hf = myFigure(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

figPos=[1186  343   680   556]; %[1432 421 1000 750];
col='w';

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case{'pos','position','figpos'}
            i=i+1;
            figPos=varargin{i};
            case{'color','col','background'}
            i=i+1;
            figPos=varargin{i};
        otherwise
            error('%s - Illegal argument: %s',current_function().upper,varargin{i});
    end
    i=i+1;
end



hf=figure;

set(hf,'position',figPos,'Color',col)


end

