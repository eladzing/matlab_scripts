function ds = ds_sphere(MPc)
% Returns a 256.^3 cube with the differential area elements (dA) matching
% a spherical cube of MPc size (See: cart2sphere).
%
% @param MPc  The required Mpc resolution (1, 2, 4, or 8)
%
% @returns     The differential area cube
%
narginchk(1,1);

[mesh_r mesh_phi mesh_theta] = sphere_grid(MPc);

n = 256*256;
ds = (mesh_r.^2)*4*pi/n;

