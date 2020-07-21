function res = nfwA(x,c)
%aNFW  NFW Auxilary function 
%   Auxilary function used for finding the Mass profile for a NFW halo 
%   also used for normalization 
%   x = r/Rvir  is in units of Rvir.
%   c is the concentration parameter;

xx=c.*x;

res=log(1+xx)-xx./(1+xx);

end

