global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

MPSec = 1;
h = 0.7;

ds = ds_sphere(MPSec);
RRHOG = RHOG_sphere(MPSec);
RRHOTOT = RHOTOT_sphere(MPSec);
TT = T_sphere(MPSec);

RHOG_Profile   = squeeze(sum(sum(RRHOG,3),2))'/(256^2);
RHOTOT_Profile = squeeze(sum(sum(RRHOTOT,3),2))'/(256^2);
T_Profile = squeeze(sum(sum(TT.*RRHOG,3),2))'./squeeze(sum(sum(RRHOG,3),2))';

[mesh_r mesh_phi mesh_theta] = sphere_grid(MPSec);
R_Profile = squeeze(mesh_r(:,1,1))';

%we would like to integrate the Rhos to we multiply by dr
dr = R_Profile(2)-R_Profile(1);

%%% WARNING: We miss the inner sphere R<5, profile is shifted by constant
MG_Profile   = cumsum(dr*squeeze(sum(sum(ds.*RRHOG,3),2))');
MTOT_Profile = cumsum(dr*squeeze(sum(sum(ds.*RRHOTOT,3),2))');
S_Profile = T_Profile ./ ((RHOG_Profile.^(2/3)));

% Plots!
%figure(1); semilogy(R_Profile, T_Profile, '*-'); title('T Profile'); xlabel('R (MPc)'); ylabel('log(T(R))'); saveas(gcf, getresultsdir(sprintf('T Profile (%d MPc).png',MPSec)));
%figure(2); plot(R_Profile, S_Profile, '*-'); title('Entropy Profile'); xlabel('R (MPc)'); ylabel('Entropy'); saveas(gcf, getresultsdir(sprintf('S Profile (%d MPc).png',MPSec)));
%figure(3); plot(R_Profile, MG_Profile, '*-'); title('Total Gas Mass Profile'); xlabel('R (MPc)'); ylabel('Total Gas Mass (Solar Masses)'); saveas(gcf, getresultsdir(sprintf('Cumulative Gas Mass Profile (%d MPc).png',MPSec)));
%figure(4); plot(R_Profile, MTOT_Profile, '*-'); title('Total Mass Profile'); xlabel('R (MPc)'); ylabel('Total Mass (Solar Masses)'); saveas(gcf, getresultsdir(sprintf('Cumulative Total Mass Profile (%d MPc).png',MPSec)));
% TODO

%figure(5); loglog(R_Profile, T_Profile, '*-'); title('T Profile'); xlabel('log(R) (MPc)'); ylabel('log(T(R))'); saveas(gcf, getresultsdir(sprintf('T Profile (%d MPc) loglog.png',MPSec)));
%figure(6); loglog(R_Profile, S_Profile, '*-'); title('Entropy Profile'); xlabel('log(R) (MPc)'); ylabel('Entropy'); saveas(gcf, getresultsdir(sprintf('S Profile (%d MPc) loglog.png',MPSec)));
%figure(7); loglog(R_Profile, MG_Profile, '*-'); title('Total Gas Mass Profile'); xlabel('log(R) (MPc)'); ylabel('Total Gas Mass (Solar Masses)'); saveas(gcf, getresultsdir(sprintf('Cumulative Gas Mass Profile (%d MPc) loglog.png',MPSec)));
%figure(8); loglog(R_Profile, MTOT_Profile, '*-'); title('Total Mass Profile'); xlabel('log(R) (MPc)'); ylabel('Total Mass (Solar Masses)'); saveas(gcf, getresultsdir(sprintf('Cumulative Total Mass Profile (%d MPc) loglog.png',MPSec)));


%% load units
units

%%
jcos = (2.76e-30/1e3*((1e2)^3))/Msun*(MPc^3)
Average_RHOTOT_Profile = (MTOT_Profile ./ ((4*pi/3)*(R_Profile.^3)));
figure(7); loglog(R_Profile, Average_RHOTOT_Profile/jcos); title('Measuring The Virial Radius (Assuming j_{cos} = 2.76\cdot10^{-30} g/cc)'); xlabel('R (MPc)'); ylabel('Average Density within R / j_{cos}');

%%
%R_Profile(find(Average_RHOTOT_Profile./jcos < 200,1))
%We take ridx or ridx-1 whichever is closer to 200

RVIR_IDX = {};
RVIR_IDX{4} = 197;
RVIR_IDX{8} = 98;

RVIR = R_Profile(RVIR_IDX{MPSec}) %Mpc
MVIR = MTOT_Profile(RVIR_IDX{MPSec}) %Msun
VVIR = sqrt(G*(MVIR*Msun)/(RVIR*MPc))*1e-3 %km/sec
TVIR = ((VVIR*1e3)^2)*m/2/k %kelvin
RHOVIR = MVIR/((4*pi/3)*(RVIR^3));
% save('virial4', 'RVIR_IDX', 'RVIR', 'MVIR', 'VVIR', 'TVIR', 'RHOVIR');

%figure(9);  loglog(R_Profile/RVIR, T_Profile/TVIR, '*-'); title('T Profile'); xlabel('r / R_{v}'); ylabel('T / T_{v}'); saveas(gcf, getresultsdir(sprintf('T Profile (%d MPc) loglog.png',MPSec)));
%figure(10); loglog(R_Profile/RVIR, S_Profile, '*-'); title('Entropy Profile'); xlabel('r / R_{v}'); ylabel('Entropy'); saveas(gcf, getresultsdir(sprintf('S Profile (%d MPc) loglog.png',MPSec)));
%figure(11); loglog(R_Profile/RVIR, MG_Profile/MG_Profile(RVIR_IDX{MPSec}), '*-'); title('Total Gas Mass Profile'); xlabel('r / R_{v}'); ylabel('Mg / Mg_{v}'); saveas(gcf, getresultsdir(sprintf('Cumulative Gas Mass Profile (%d MPc) loglog.png',MPSec)));
%figure(12); loglog(R_Profile/RVIR, MTOT_Profile/MVIR, '*-'); title('Total Mass Profile'); xlabel('r / R_{v}'); ylabel('M / M_{v}'); saveas(gcf, getresultsdir(sprintf('Cumulative Total Mass Profile (%d MPc) loglog.png',MPSec)));
%figure(13); loglog(R_Profile/RVIR, RHOG_Profile/RHOVIR, '*-'); title('Gas Denstiy Profile'); xlabel('r / R_{v}'); ylabel('\rho_{g} / \rho_{v}'); %saveas(gcf, getresultsdir(sprintf('Cumulative Gas Mass Profile (%d MPc) loglog.png',MPSec)));
%figure(14); loglog(R_Profile/RVIR, RHOTOT_Profile/RHOVIR, '*-'); title('Total Density Profile'); xlabel('r / R_{v}'); ylabel('\rho / \rho_{v}'); %saveas(gcf, getresultsdir(sprintf('Cumulative Total Mass Profile (%d MPc) loglog.png',MPSec)));

Normalized_S_Profile = (T_Profile/TVIR) ./ (((RHOG_Profile/RHOVIR).^(2/3)));
figure(15); loglog(R_Profile'/RVIR, [T_Profile'/TVIR RHOG_Profile'/RHOVIR RHOTOT_Profile'/RHOVIR Normalized_S_Profile'], '.-'); title('Profiles'); xlabel('r / R_{v}'); legend('T / T_{v}', '\rho_{g} / \rho_{v}', '\rho / \rho_{v}', 'Gas Entropy K= (T/T_{v})/(\rho/\rho_{v})^{2/3})'); xlim([min(R_Profile/RVIR) max(R_Profile/RVIR)]); saveas(gcf, getresultsdir(sprintf('Profiles (%d MPc).png',MPSec)));


%%
LEFT  = (-1) * G * (MTOT_Profile*Msun) ./ ((R_Profile*MPc).^2);
LEFT = (LEFT(2:end) + LEFT(1:end-1))/2;
Fixed_RHOG_Profile = RHOG_Profile * Msun / ((MPc)^3);
RIGHT_BRACKETS = (k / m) * (Fixed_RHOG_Profile .* T_Profile);
Fixed_RHOG_Profile = (Fixed_RHOG_Profile(2:end) + Fixed_RHOG_Profile(1:end-1))/2;
RIGHT = (1./Fixed_RHOG_Profile) .* ((RIGHT_BRACKETS(2:end)-RIGHT_BRACKETS(1:end-1)) / (dr*MPc));

%%
R2_Profile = (R_Profile(2:end) + R_Profile(1:end-1))/2;
%figure; plot(R2_Profile/RVIR, abs(RIGHT./LEFT), '.-'); line([min(R2_Profile/RVIR) max(R2_Profile/RVIR)], [1 1]); title('Hydrostatic Equilibrium'); ylabel('Pressure Term / Gravitational Term'); xlabel('r / R_{v}'); saveas(gcf, getresultsdir(sprintf('Hydrostatic Radio (%d MPSec).png',MPSec)));
%figure; loglog(R2_Profile/RVIR, abs([RIGHT' LEFT']), '.-'); title('Hydrostatic Equilibrium'); xlabel('r / R_{v}'); ylim([0,1]*1e-9); legend('Pressure Term', 'Gravitational Term'); saveas(gcf, getresultsdir(sprintf('Hydrostatic (%d MPSec).png',MPSec)));
figure; plot(R2_Profile/RVIR, abs(RIGHT./LEFT), '.-'); line([min(R2_Profile/RVIR) max(R2_Profile/RVIR)], [1 1]); title('Hydrostatic Equilibrium'); ylabel('Pressure Term / Gravitational Term'); xlabel('r / R_{v}'); saveas(gcf, getresultsdir(sprintf('Hydrostatic Radio (%d MPSec).png',MPSec)));
figure; loglog(R2_Profile/RVIR, abs([RIGHT' LEFT']), '.-'); title('Hydrostatic Equilibrium'); xlabel('r / R_{v}'); legend('Pressure Term', 'Gravitational Term'); saveas(gcf, getresultsdir(sprintf('Hydrostatic (%d MPSec).png',MPSec)));

%% (ONLY WORKS AFTER RUNNING work_flux_sphere
figure(22); 
plot(R2_Profile/RVIR, abs(RIGHT./LEFT), '.-', RR/RVIR, [abs(Mout_DOT./Min_DOT)]', 'r.-');
legend('Pressure Term / Gravitational Term', 'Mdot_{out} / Mdot_{in}');
line([min(R2_Profile/RVIR) max(R2_Profile/RVIR)], [1 1]);
xlabel('r / R_{v}');
title('Hydrostatic Equilibrium vs. Flux');
ylim([0, 1.6]);
saveas(gcf, getresultsdir(sprintf('Hydrostatic Radio vs. Flux (%d MPSec).png',MPSec)));

