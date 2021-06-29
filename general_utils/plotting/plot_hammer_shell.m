function plot_hammer_shell(shell,theta,phi,varargin)
%% function for plotting a spherical shell on a hammar map
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



%meshgrid of the projection, 360 deg horizontaly, 180 vertically.
theta2=(theta+pi).*180/pi;
phi2=(phi+0.5.*pi).*180/pi;

dth=360/nc;
dph=180/(0.5*nc);
[thetaI,phiI] = meshgrid(dth:dth:360,180:-dph:dph);
%[thetaI,phiI] = meshgrid(1:360,1:180);

% interpolate the data at the positions given by thetaI and phiI
dataitoff = double(interp2(theta2,phi2,shell,thetaI,phiI));

%command to draw the hammer projection axis
axesm('hammer','Grid', 'on');

% plot the dataitoff on the projection axis. 
meshm(dataitoff,[nc/360,90,180]);
if barflag
    colorbar;
end


