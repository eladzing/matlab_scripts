function [dydx,xNew]=deriveLog1(y,x)
%% simple function for numerical evaluation of derivative 
% the function evaluates the derivative at indices 2:end-1
% and returns the new x vector

dx=diff(log10(x));
dy=diff(log10(y));

dydx=(dy(2:end)+dy(1:end-1))./(dx(1:end-1)+dx(2:end));

xNew=x(2:end-1);