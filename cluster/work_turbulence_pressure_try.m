new_env(3)
MPSec = 8;
[vx, vy, vz, v1,v2,v3] = V_Vcm_new(MPSec);

%[bins_x_vx bins_vx] = calc_kolmogorov(vx);
%[bins_x_vy bins_vy] = calc_kolmogorov(vy);
%[bins_x_vz bins_vz] = calc_kolmogorov(vz);

%%
%%loglog(bins_x_vx, bins_vx,bins_x_vy, bins_vy,bins_x_vz, bins_vz);
%loglog(1:length(bins_vx), [bins_vx;bins_vy;bins_vz]);
%%loglog(bins_x_vx, [bins_vx;bins_vy;bins_vz], '-', 10:80, (10:80).^(-2.5).*(10^(9.4)),'--');

%%
%we see a cutoff at k=9. we use it to build a mask for the fourier
%transform:
small_len = length(vx);
grid_vec = (1:small_len)-((small_len+1)/2);
[xxx,yyy,zzz] = meshgrid(grid_vec,grid_vec,grid_vec);
radii = sqrt(xxx.^2+yyy.^2+zzz.^2);
clear xxx yyy zzz
%mask = bitand(radii > 9, radii < 81);
mask = (radii > 9);
clear radii

%%
vx2 = ifftn(ifftshift(fftshift(fftn(vx)).*double(mask)));
vy2 = ifftn(ifftshift(fftshift(fftn(vy)).*double(mask)));
vz2 = ifftn(ifftshift(fftshift(fftn(vz)).*double(mask)));

%% sanity check
%[bins_x_vx bins_vx] = calc_kolmogorov(vx2);
%[bins_x_vy bins_vy] = calc_kolmogorov(vy2);
%[bins_x_vz bins_vz] = calc_kolmogorov(vz2);
%loglog(bins_x_vx, [bins_vx;bins_vy;bins_vz], '-', 10:80, (10:80).^(-2.5).*(10^(9.4)),'--');

%%
v_sq = vx.^2 + vy.^2 + vz.^2;
new_v_sq = abs(vx2).^2 + abs(vy2).^2 + abs(vz2).^2;
%%
units
%this line calculates the pressure locally. instead we can calculate only
%the V_sq profile and multiply it by RHOG profile
%pressure = 0.5*new_v_sq.*1e3.*RHOG(1)*(Ms/((Mpc)^3)); %in MKS!
pressure = 0.5*new_v_sq.*1e6;
pressure_sphere = cart2sphere(pressure);
pressure_profile  = squeeze(sum(sum(pressure_sphere,3),2))'/(256^2);

%%
% save('turbulence.mat', 'pressure_sphere', 'pressure_profile', 'pressure', 'v_sq', 'new_v_sq');
%% Run this after running work_profiles_sphere_new
%%load turbulence
%Fixed_RHOG_Profile = RHOG_Profile * Ms / ((Mpc)^3);%
%figure(1);semilogy(R_Profile/RVIR,[RIGHT_BRACKETS;pressure_profile.*Fixed_RHOG_Profile]);xlabel('r/R_v');ylabel('Pressure (N/m^2)');title('Comparing pressure terms'); legend('Thermal Pressure = k/m*\rho_g(r)*T(r)', 'Turbulent Pressure = 1/2*\rho_g(r)*\sigma^2');saveas(gcf,getresultsdir('pressure terms.png'));

%% Read profiles
global NCELL 
global hub
R_Profile=MPSec/NCELL/2:MPSec/NCELL/2:MPSec/2;
R_Profile=R_Profile./hub;

RVIR=get_rvir;

MTOT_Profile = read_MTOT_Profile(R_Profile);
[RHOG_Profile,~] = read_RHO_Profiles(R_Profile);
T_Profile = read_T_Profile(R_Profile);
dr= max(diff(R_Profile));

%%

units

LEFT  = (-1) * G * (MTOT_Profile*Ms) ./ ((R_Profile*Mpc).^2);
LEFT = (LEFT(2:end) + LEFT(1:end-1))/2;
Fixed_RHOG_Profile = RHOG_Profile * Ms / ((Mpc)^3);

RIGHT_BRACKETS = (kb / mp) * (Fixed_RHOG_Profile .* T_Profile);
NEW_RIGHT_BRACKETS = ((kb / mp) * (Fixed_RHOG_Profile .* T_Profile)) + (pressure_profile.*Fixed_RHOG_Profile);

Fixed_RHOG_Profile = (Fixed_RHOG_Profile(2:end) + Fixed_RHOG_Profile(1:end-1))/2;

RIGHT = (1./Fixed_RHOG_Profile) .* ((RIGHT_BRACKETS(2:end)-RIGHT_BRACKETS(1:end-1)) / (dr*Mpc));
NEW_RIGHT = (1./Fixed_RHOG_Profile) .* ((NEW_RIGHT_BRACKETS(2:end)-NEW_RIGHT_BRACKETS(1:end-1)) / (dr*Mpc));



%%
R2_Profile = (R_Profile(2:end) + R_Profile(1:end-1))/2;
figure; loglog(R2_Profile/RVIR, [abs(RIGHT./LEFT)' abs(NEW_RIGHT./LEFT)'], '.-'); line([min(R2_Profile/RVIR) max(R2_Profile/RVIR)], [1 1]); title('Hydrostatic Equilibrium inc. Turbulent Pressure'); ylabel('Pressure Term / Gravitational Term'); xlabel('r / R_{v}'); legend('not including turbulentce', 'including turbulence' );%saveas(gcf, getresultsdir(sprintf('Hydrostatic Ratio with Turbulence (%d MPSec).png',MPSec)));
figure; loglog(R2_Profile/RVIR, abs([RIGHT' NEW_RIGHT' LEFT']), '.-'); title('Hydrostatic Equilibrium inc. Turbulent Pressure'); xlabel('r / R_{v}'); legend('Pressure Term (not inc. Turbulence)', 'Pressure Term (inc. Turbulence)', 'Gravitational Term'); %saveas(gcf, getresultsdir(sprintf('Hydrostatic with Turbulence (%d MPSec).png',MPSec)));
