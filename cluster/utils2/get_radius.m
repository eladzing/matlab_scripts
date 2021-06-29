function [rad, radR] = get_radius(figHandle,np)
%GET_RADIUS - find radius of point on figure, also in units of Rvir
%   get_radius returns the radius in the plain of the picture after user
%   interacitvely chooses points on figure figHandle.
%   The radius is returned in physical units and in units of the virial
%   radius. The plots must originally be in proper units (not comoving) to
%   allow the comparision to Rvir.

figure(figHandle)

if exist('np','var')
    [xx, yy]=ginput(np);
else
    [xx, yy]=ginput();
end

global hub

xx=xx./hub;
yy=yy./hub;

rad=sqrt(xx.^2+yy.^2);
radR=rad./get_rvir;

end

