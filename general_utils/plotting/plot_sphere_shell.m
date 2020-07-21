function plot_sphere_shell(shell,theta,phi,varargin)
%% function for plotting a spherical shell on a spherical map
% shell is a 2d array of values extracted from one of the spherical cubes
% theta - array of theta values (actually phi in normal convention)
% phi - array of phi values (actually theta in normal convention

nc=360;
barflag=false;

i=1;
while i<=nargin-3
    switch varargin{i}
        case 'nc'
            i=i+1;
            nc=varargin{i};
        case {'colorbar','bar'}
            barflag=true;
        otherwise
            error('plot_hammer_shell: illegal argument');
    end
    i=i+1;
end

theta2=(theta+pi);
phi2=(phi+0.5.*pi);

%meshgrid of the projection, 360 deg horizontaly, 180 vertically.

dth=2*pi/nc;
dph=pi/(0.5*nc);
[thetaI,phiI] = meshgrid(dth:dth:2*pi,pi:-dph:dph);
% %[thetaI,phiI] = meshgrid(1:360,1:180);
%
% interpolate the data at the positions given by thetaI and phiI
ishell = double(interp2(theta2,phi2,shell,thetaI,phiI));

zz=cos(phiI);
yy=sin(thetaI).*sin(phiI);
xx=cos(thetaI).*sin(phiI);

surf(xx,yy,zz,double(ishell),'edgecolor','none')
if barflag
    colorbar;
end


