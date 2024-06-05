function [mesh_r, mesh_phi, mesh_theta] = sphere_grid(boxSize,NCELL)
% Utility function to generate the mesh grids for the spherical coordinates returned by cart2sphere().
% Coordinates are (R,Phi,Theta) such that:
%    R is in units of Mpc (between 0.5*(mpc/h)/256 and 128*(mpc/h)/256)
%    -Pi/2 <= Phi   <= Pi/2 (in non-linear spacing to impose a constant ds
%                          per radius)
%    -Pi   <= Theta <= Pi
%
% @param mpc  The required Mpc resolution (1, 2, 4, or 8)
% @returns    The three mesh grids for the R, Phi, and Theta coordinates.
%
% argument boxx is the size of the box in Mpc/h

%narginchk(1,1);

n = NCELL-1;
theta = (-n:2:n)/n*pi;
phi = acos((-n:2:n)'/n)-pi/2;






[mesh_phi, mesh_r, mesh_theta] = meshgrid(phi, (0.5:0.5:NCELL/2).*boxSize./NCELL, theta);