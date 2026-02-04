function [dydx,xNew]=derive1(y,x,varargin)
%% simple function for numerical evaluation of derivative 
% the function evaluates the derivative at indices 2:end-1
% and returns the new x vector
% user can specify along which dimension to do the derivative (up to 2) 

sx=size(x);sy=size(y);

if length(sx)>2 || length(sy) >>2
    error('%s - function not equipped to handle arrays larger than 2D \n',current_function().upper)
end

if sx~=sy
    error('%s - input array do not match in size: size y= %i %i, size x = %i, %i',...
        current_function().upper,sy(1),sy(2),sx(1),sx(2));
end

ddim=[];

%% parse argument 
i=1;
while(i<=length(varargin))
    switch(lower(varargin{i}))
        case{'dim','dimension'}
            i=i+1;
            ddim=varargin{i};

            if ddim~=1 && ddim~=2
                error('%s - dimension should be 1 or 2, not: %i',current_function().upper,varargin{i})
            end

        otherwise   
                error('%s - Illegal argument: %s',current_function().upper,varargin{i})
    end
    i=i+1;
end

%% define dimension to derive along

if isempty(ddim) || (any(sx==1) && any(sy==1))
    if ~any(sx==1) || ~any(sy==1)
        error('%s - must declare dimension if array is not 1D',current_function().upper);
    end

    dx=diff(x);
    dy=diff(y);

    dydx=(dy(2:end)+dy(1:end-1))./(dx(1:end-1)+dx(2:end));

    xNew=x(2:end-1);

else

    dx=diff(x,1,ddim);
    dy=diff(y,1,ddim);

    if ddim==1
        dydx=(dy(2:end,:)+dy(1:end-1,:))./(dx(1:end-1,:)+dx(2:end,:));
        xNew=x(2:end-1,:);

    elseif ddim==2
        dydx=(dy(:,2:end)+dy(:,1:end-1))./(dx(:,1:end-1)+dx(:,2:end));
        xNew=x(:,2:end-1);
    else
        error('%s - We should not be able to get here...',current_function().upper)
    end
    

