function res=bfunc_derivative(r,nu)
%% auxilary function to calculate the 'bessel' part of the  
% gravitational acceleration of an exponential disk.
% r is in units of the scale radius rd 
% nu is the ratio between rd and rg (star & gas scale radii)

if ~exist('nu','var')
    nu=1;
end

i1=besseli(1,0.5.*r.*nu);
i0=besseli(0,0.5.*r.*nu);
k0=besselk(0,0.5.*r.*nu);
k1=besselk(1,0.5.*r.*nu);

%res=besseli(0,0.5.*r*nu).*besselk(0,0.5.*r.*nu)-besseli(1,0.5.*r.*nu).*besselk(1,0.5.*r.*nu);

res=nu.*(i1.*k0-i0.*k1)+2.*i1.*k1./r;

%res=r.*nu.*(besseli(0,0.5.*r*nu).*besselk(0,0.5.*r.*nu)-besseli(1,0.5.*r.*nu).*besselk(1,0.5.*r.*nu));

