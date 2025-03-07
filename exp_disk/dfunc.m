function res=dfunc(r,nu)
%% auxilary function to calculate the 'bessel' part of the  
% gravitational potential of an exponential disk.
% r is in units of the scale radius rd 
% nu is the ratio between rd and rg (star & gas scale radii)

if ~exist('nu','var')
    nu=1;
end


res=r.*(besseli(0,0.5.*r*nu).*besselk(1,0.5.*r.*nu)-besseli(1,0.5.*r.*nu).*besselk(0,0.5.*r.*nu));

