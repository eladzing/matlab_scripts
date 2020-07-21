%by hand I calculated the entropy threshold bi-modality in different radii
%This line is taken from work_bimodality
S8 = load('../data/S_sphere_8','res'); S8 = S8.res;
ds8 = ds_sphere(8);
RG8 = load('../data/RHOG_sphere_8','res'); RG8 = RG8.res;
place = [110:150];[range vals] = hist_density((S8(place,:,:)),RG8(place,:,:).*ds8(place,:,:),100,0,0.35); bar(range,vals);
%
%where place was varying between 20:30 until 100:110. ranges below or above
%where too noisy.
%we come with the following estimation
%[mesh_r mesh_phi mesh_theta] = sphere_grid(8);
%xx in Mpc and yy in the entropy values
xx = [25,35,45,55,65,75,85,95,105]/256*4/0.7; yy = [0.05,0.07,0.09,0.1,0.11,0.13,0.145,0.17,0.19];

