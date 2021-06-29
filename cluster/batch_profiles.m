function found_rvir = batch_profiles(MPSec, halopath, clustername, load_virial)

rhog_sp = RHOG_sphere(MPSec);
rhotot_sp = RHOTOT_sphere(MPSec);
t_sp = T_sphere(MPSec);
ds = ds_sphere(MPSec);

RHOG_Profile   = double(squeeze(sum(sum(rhog_sp,3),2))'/(256^2));
RHOTOT_Profile = double(squeeze(sum(sum(rhotot_sp,3),2))'/(256^2));
T_Profile = double(squeeze(sum(sum(t_sp.*rhog_sp,3),2))'./squeeze(sum(sum(rhog_sp,3),2))');

[mesh_r mesh_phi mesh_theta] = sphere_grid(MPSec);
clear mesh_phi mesh_theta
R_Profile = double(squeeze(mesh_r(:,1,1))');

%we would like to integrate the Rhos to we multiply by dr
dr = double(R_Profile(2)-R_Profile(1));

%%% WARNING: We miss the inner sphere R<5, profile is shifted by constant
%MG_Profile   = cumsum(dr*squeeze(sum(sum(ds.*rhog_sp,3),2))');
MTOT_Profile = cumsum(dr*squeeze(sum(sum(double(ds.*rhotot_sp),3),2))');
S_Profile = T_Profile ./ ((RHOG_Profile.^(2/3)));


%% load units
units

%%
jcos = (2.76e-30/1e3*((1e2)^3))/Msun*(MPc^3);
Average_RHOTOT_Profile = (MTOT_Profile ./ ((4*pi/3)*(R_Profile.^3)));

%% Find virial quantities

found_rvir = 1;

if (load_virial ~= 1)
    Average_RHOTOT_Profile = Average_RHOTOT_Profile./jcos;
    RVIR_IDX = find(Average_RHOTOT_Profile < 200,1);
    if (length(RVIR_IDX) == 0)
        found_rvir = 0;
        return
    end
    
    %We take ridx or ridx-1 whichever is closer to 200
    if (abs(Average_RHOTOT_Profile(RVIR_IDX)-200) > abs(Average_RHOTOT_Profile(RVIR_IDX-1)-200))
        RVIR_IDX = RVIR_IDX-1;
    end

    RVIR = R_Profile(RVIR_IDX)

    MVIR = MTOT_Profile(RVIR_IDX); %Msun
    VVIR = sqrt(G*(MVIR*Msun)/(RVIR*MPc))*1e-3; %km/sec
    TVIR = ((VVIR*1e3)^2)*m/2/k; %kelvin
    RHOVIR = MVIR/((4*pi/3)*(RVIR^3));
else
    load(sprintf('%s/virial%d.mat',halopath, 8), 'RVIR', 'TVIR', 'RHOVIR', 'VVIR', 'MVIR');
end

Normalized_S_Profile = (T_Profile/TVIR) ./ (((RHOG_Profile/RHOVIR).^(2/3)));
figure; loglog(R_Profile'/RVIR, [T_Profile'/TVIR RHOG_Profile'/RHOVIR RHOTOT_Profile'/RHOVIR Normalized_S_Profile'], '.-'); 
title('Profiles'); xlabel('r / R_{v}'); legend('T / T_{v}', '\rho_{g} / \rho_{v}', '\rho / \rho_{v}', 'Gas Entropy K= (T/T_{v})/(\rho/\rho_{v})^{2/3})'); 
xlim([min(R_Profile/RVIR) max(R_Profile/RVIR)]); saveas(gcf, getresultsdir(sprintf('%s - Profiles (%d MPc).png',clustername, MPSec)));

%% Calculate HSE
LEFT  = (-1) * G * (MTOT_Profile*Msun) ./ ((R_Profile*MPc).^2);
LEFT = (LEFT(2:end) + LEFT(1:end-1))/2;
Fixed_RHOG_Profile = RHOG_Profile * Msun / ((MPc)^3);
RIGHT_BRACKETS = (k / m) * (Fixed_RHOG_Profile .* T_Profile);
Fixed_RHOG_Profile = (Fixed_RHOG_Profile(2:end) + Fixed_RHOG_Profile(1:end-1))/2;
RIGHT = (1./Fixed_RHOG_Profile) .* ((RIGHT_BRACKETS(2:end)-RIGHT_BRACKETS(1:end-1)) / (dr*MPc));

R2_Profile = (R_Profile(2:end) + R_Profile(1:end-1))/2;
figure; plot(R2_Profile/RVIR, abs(RIGHT./LEFT), '.-'); line([min(R2_Profile/RVIR) max(R2_Profile/RVIR)], [1 1]); title('Hydrostatic Equilibrium'); ylabel('Pressure Term / Gravitational Term'); xlabel('r / R_{v}'); saveas(gcf, getresultsdir(sprintf('%s - Hydrostatic Ratio (%d MPSec).png',clustername, MPSec)));
figure; loglog(R2_Profile/RVIR, abs([RIGHT' LEFT']), '.-'); title('Hydrostatic Equilibrium'); xlabel('r / R_{v}'); legend('Pressure Term', 'Gravitational Term'); saveas(gcf, getresultsdir(sprintf('%s - Hydrostatic (%d MPSec).png', clustername, MPSec)));

%% Save Data
virial_filename = sprintf('%s/virial%d.mat',halopath, MPSec);
if (load_virial ~= 1)
    save(virial_filename, 'MVIR', 'RVIR', 'VVIR', 'RVIR_IDX', 'RHOVIR', 'TVIR', 'MTOT_Profile', 'T_Profile', 'RHOG_Profile','RHOTOT_Profile', 'R_Profile', 'Average_RHOTOT_Profile' , 'Normalized_S_Profile', 'S_Profile', 'RIGHT', 'LEFT', 'R2_Profile');
else
    save(virial_filename, 'MTOT_Profile', 'T_Profile', 'RHOG_Profile','RHOTOT_Profile', 'R_Profile', 'Average_RHOTOT_Profile' , 'Normalized_S_Profile', 'S_Profile', 'RIGHT', 'LEFT', 'R2_Profile');
end