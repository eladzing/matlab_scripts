function [rad, radR] = get_radius3(hf1,hf2,hf3)
%GET_RADIUS - find the full 3D radius of a point based on 3 projections 
%   get_radius3 returns the full radius based on 3 figures of the 3 projections 
%   interacitvely chooses points on figure figHandle. 
%   The radius is returned in physical units and in units of the virial
%   radius. The plots must originally be in proper units (not comoving) to
%   allow the comparision to Rvir. 


[rd(1),~] = get_radius(hf1,1);
[rd(2),~] = get_radius(hf2,1);
[rd(3),~] = get_radius(hf3,1);

rad=sqrt(0.5.*sum(rd.^2));

radR=rad/get_rvir;

end

