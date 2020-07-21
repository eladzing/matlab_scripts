function hub_sp = hubble_sphere(MPc) 
%%% returns a spherical cube containing the hubble flow 
%%% velocity in km/sec from the center of the box 

[mesh_r mesh_phi mesh_theta] = sphere_grid(MPc);
clear mesh_phi mesh_theta;

global hub

hub_sp=hub.*100.*mesh_r;
