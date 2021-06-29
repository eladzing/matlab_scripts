env

MPSec = 8;

[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);
[Vrr cosa] = Vr_mesh(MPSec);

mask = (cosa<0);
clear cosa Vrr VcmX VcmY VcmZ

Vxx(~mask) = 0;
Vyy(~mask) = 0;
Vzz(~mask) = 0;

clear mask 

%%
% DILUTE = 3;
% 
% [startY1 startX1 startZ1] = meshgrid([1:DILUTE:256], [1 256], [1:DILUTE:256]);
% [startY2 startX2 startZ2] = meshgrid([1 256], [1:DILUTE:256], [1:DILUTE:256]);
% [startY3 startX3 startZ3] = meshgrid([1:DILUTE:256], [1:DILUTE:256], [1 256]);
% 
% [meshY, meshX, meshZ] = meshgrid(1:size(Vyy,1), 1:size(Vxx,2), 1:size(Vzz,3));
% 
% verts = stream3(meshY, meshX, meshZ, Vyy, Vxx, Vzz, [startY1(:) startY2(:) startY3(:)], [startX1(:) startX2(:) startX3(:)], [startZ1(:) startZ2(:) startZ3(:)]);

%%

[meshY, meshX, meshZ] = meshgrid(1:size(Vyy,1), 1:size(Vxx,2), 1:size(Vzz,3));
SPHERE_RES = 256;
[Sx Sy Sz] = uni_sphere(SPHERE_RES);
Sx = (Sx*127.5)+128.5;
Sy = (Sy*127.5)+128.5;
Sz = (Sz*127.5)+128.5;
verts = stream3(meshY, meshX, meshZ, Vyy, Vxx, Vzz, Sy, Sx, Sz);

% n = SPHERE_RES^2;
% ds = (mesh_r.^2)*4*pi/n;


%mask = f(verts, 1.5/8*256);
%%
grid_vec = (1:256)-((256+1)/2);
[xxx,yyy,zzz] = meshgrid(grid_vec,grid_vec,grid_vec);
radii_square = (xxx.^2+yyy.^2+zzz.^2);
clear xxx yyy zzz grid_vec
h=0.7;
radii_square = radii_square/(256^2)*((8/h)^2);
%%
%TODO: select only streamlines that did not pass a shock!

clear Vxx Vyy Vzz Sx Sy Sz meshX meshY meshZ

%disregard streamlines that don't get as close as Rv
bins_factor = 50
[f_results f_stds f_min f_max f_sum f_sum_count] = stream_profile_bins(verts, 1.5/8*256, RHOG(8).*Vr(8).*radii_square*4*pi/(SPHERE_RES^2), length(verts), bins_factor);

%%
figure;
load profiles8
XX = [1:length(f_results)]/bins_factor;
f_sum(f_sum==0)=NaN;
plot(XX, f_sum/abs(used_mdot),'.-')
xlabel('r/R_{v}')
ylabel('flux')
legend('flux of stream lines', 'STD bars','STD bars', 'Location', 'NorthWest');
title('Flux of streamlines vs. cube')
saveas(gcf, getresultsdir('Flux of streamlines vs. cube (with stds).png'))

%% after running work_flux_sphere
plot(XX, f_sum/abs(used_mdot),'.-', RR/RVIR, Mcold_DOT/abs(used_mdot))
xlabel('r/R_{v}')
ylabel('flux Mdot/Mdot_v')
legend('stream lines', 'cold flows', 'Location', 'NorthWest');
title('Flux of streamlines vs. cold flows')
saveas(gcf, getresultsdir('Flux of streamlines vs. cold flows.png'))
