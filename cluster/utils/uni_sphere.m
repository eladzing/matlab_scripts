function [xx,yy,zz] = uni_sphere(varargin)
%UNI_SPHERE Generate sphere.
% IMPORTANT: THIS IS A MODIFIED VERSION OF MATLAB'S SPHERE() FUNCTION
%            We modified the Phi coordinate spacing to enforce a constant
%            ds per radius, and generate a UNIform sphere.
%
%   [X,Y,Z] = UNI_SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURF(X,Y,Z) produces a unit sphere.
%
%   [X,Y,Z] = UNI_SPHERE uses N = 20.
%
%   UNI_SPHERE(N) and just SPHERE graph the sphere as a SURFACE
%   and do not return anything.
%
%   UNI_SPHERE(AX,...) plots into AX instead of GCA.
%
%   See also ELLIPSOID, CYLINDER.

%   Clay M. Thompson 4-24-91, CBM 8-21-92.
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.8.4.1 $  $Date: 2002/09/26 01:55:25 $

% Parse possible Axes input
narginchk(0,2);
[cax,args,nargs] = axescheck(varargin{:});

n = 20;
if nargs > 0, n = args{1}; end

% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = acos((-n:2:n)'/n)-pi/2; %%%%%%(-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;

x = cosphi*cos(theta);
y = cosphi*sintheta;
z = sin(phi)*ones(1,n+1);

if nargout == 0
    cax = newplot(cax);
    surf(x,y,z,'parent',cax)
else
    xx = x; yy = y; zz = z;
end
