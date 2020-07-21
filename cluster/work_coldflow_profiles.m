env

%%
load ('/home/alf/giorae/results/2008-01-10/coldflows/parametric 2008-01-10 16-20-49.mat') %comment this line for interactive selection

MPSec = 8;

%%
load virial8
[mesh_r mesh_phi mesh_theta] = sphere_grid(MPSec);
R_Profile = squeeze(mesh_r(:,1,1))';
clear mesh_r mesh_phi mesh_theta

%%
rho_log = log10(RHOG_sphere(MPSec));
temp_log = log10(T_sphere(MPSec));
coldmask = inpolygon(temp_log,rho_log,xv,yv);

clear rho_log temp_log

%% density profile
RRHOG = RHOG_sphere(MPSec);
[pcube_profile masked_profile notmasked_profile] = calc_masked_profile(MPSec, RRHOG, coldmask);

figure;
loglog(R_Profile/RVIR, [pcube_profile;masked_profile;notmasked_profile]);
xlabel('r/R_{v}')
ylabel('\rho')
legend('Cube''s profile', 'Cold-flows profile','other flows profile', 'Location', 'NorthEast');
title('Density profile of cold flows')
saveas(gcf, getresultsdir('Density Profile of cold flows.png'))

%% temperature profile
TT = T_sphere(MPSec);
[pcube_profile masked_profile notmasked_profile] = calc_masked_profile(MPSec, TT, coldmask);

figure;
loglog(R_Profile/RVIR, [pcube_profile;masked_profile;notmasked_profile]);
xlabel('r/R_{v}')
ylabel('T (Kelvin)')
legend('Cube''s profile', 'Cold-flows profile','other flows profile', 'Location', 'NorthEast');
title('Temperature profile of cold flows')
saveas(gcf, getresultsdir('Temperature Profile of cold flows.png'))

%% temperature profile (wheighted by density)
[pcube_profile_w masked_profile_w notmasked_profile_w] = calc_masked_profile(MPSec, TT, coldmask, RRHOG);
figure;
loglog(R_Profile/RVIR, [pcube_profile_w;masked_profile_w;notmasked_profile_w], R_Profile/RVIR, [pcube_profile;masked_profile;notmasked_profile], '--');
xlabel('r/R_{v}')
ylabel('T (Kelvin)')
legend('Density weighted cube''s profile', 'Density weighted cold-flows profile','Density weighted other flows profile', 'Cube''s profile', 'Cold-flows profile','other flows profile', 'Location', 'SouthWest');
title('Density weighted temperature profile of cold flows')
saveas(gcf, getresultsdir('Density Weighted Temperature Profile of cold flows.png'))
