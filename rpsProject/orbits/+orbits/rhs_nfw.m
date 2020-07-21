function [xdot,ydot,vxdot,vydot] = rhs_nfw(x,y,vx,vy,host)
%RHS_NFW calcutlate the RHS of the motion equation for orbit integration
%   for an NFW halo.
%   Units assumed to be in kpc Msun km/sec 


GG=4.2997e-06;

%% set velocity as position derivative
xdot=vx;
ydot=vy;

%% set acceleration as velocity derivative

rpos=sqrt(x^2+y^2);

%acc=-1.*Units.GG.*nfw.mass(rpos)/r^3;
mm=mass(host,rpos,'kpc'); 
acc=-1.*GG.*mm/rpos^3;

vxdot=acc*x;
vydot=acc*y;

end

