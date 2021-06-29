%stack plots
global FILE_FORMAT;

base_dir = '/home/alf/eladzing/data/kravtsov/clusters';
clusternames = {'CL24csf1', 'CL6csf1', 'CL10csf1', 'CL11csf1', 'CL14csf1', ...
    'CL5csf1', 'CL7a', 'CL7csf1', 'CL9csf1', 'LCDM120_CL1_h1', 'LCDM120_CL2_h1', ...
    'LCDM120_CL3_h1', 'LCDM120_CL5_h1', 'LCDM120_CL6_h1', 'LCDM120_CL7_h1', 'LCDM120_CL8_h1', 'LCDM80_CL3_h2' };
haloformats = { '%s_a0.999L%dMpc.dat', '%s_a1.001L%dMpc.dat', '%s_a1.000L%dMpc.dat', '%s_a0.999L%dMpc.dat', ...
    '%s_a1.001L%dMpc.dat', '%s_a0.999L%dMpc.dat', '%s_a1.000L%dMpc.dat', '%s_a1.000L%dMpc.dat', '%s_a0.999L%dMpc.dat', ...
    '%s_a0.999L%dMpc.dat', '%s_a1.001L%dMpc.dat', '%s_a1.000L%dMpc.dat', '%s_a1.000L%dMpc.dat', '%s_a1.000L%dMpc.dat', ...
    '%s_a1.001L%dMpc.dat', '%s_a0.999L%dMpc.dat', '%s_a1.000L%dMpc.dat' };
h=0.7;

used_idxs = [];
for idx = 1:length(clusternames)
    halopath = sprintf('%s/%s', base_dir, clusternames{idx});
    
    try
        for MPSec = [1 2 4 8]
            virial_filename = sprintf('%s/virial%d.mat',halopath, MPSec);
            all_profiles{idx}{MPSec} = load(virial_filename);

            flux_filename = sprintf('%s/fluxes%d',halopath, MPSec);
            all_fluxes{idx}{MPSec} = load(flux_filename);

            bird_filename = sprintf('%s/bird%d.mat',halopath, MPSec);
            all_birds{idx}{MPSec} = load(bird_filename);
        end

        turb_filename = sprintf('%s/turbulence',halopath);
        all_turbs{idx} = load(turb_filename);

        used_idxs(end+1) = idx;
    catch
    end
end

units

%% plot the profiles of all clusters together
for idx = 1:length(used_idxs)
    % create entropy quiver maps
    FILE_FORMAT = sprintf('%s/%s/%s', base_dir, clusternames{used_idxs(idx)}, haloformats{used_idxs(idx)})
    texttowrite = sprintf('R_v: %d (Mpc)\nM_v: %d (M_{sun})\nT_v: %d (K)\nV_v: %d (Km/s)\n\\rho_v: %d (M_{sun}/Mpc^3)', ...
			all_profiles{used_idxs(idx)}{8}.RVIR, all_profiles{used_idxs(idx)}{8}.MVIR, all_profiles{used_idxs(idx)}{8}.TVIR, ...
			all_profiles{used_idxs(idx)}{8}.VVIR, all_profiles{used_idxs(idx)}{8}.RHOVIR);
%   work_quiver_new(clusternames{used_idxs(idx)}, 8, [], 1, texttowrite, all_profiles{used_idxs(idx)}.RVIR*h);
% 	set(gcf, 'PaperType', 'A4'); orient landscape; print -dps2 -Plwclr; print -dps2 -Plwclr


    
    % create graphs parameters
    %plot_t_params{1+2*(idx-1)} = all_profiles{used_idxs(idx)}.R_Profile'/all_profiles{used_idxs(idx)}.RVIR;
    %plot_t_params{2+2*(idx-1)} = all_profiles{used_idxs(idx)}.T_Profile'/all_profiles{used_idxs(idx)}.TVIR;
    [plot_t_params{1+2*(idx-1)} plot_t_params{2+2*(idx-1)}] = merge_multi_res_from_hash(...
        all_profiles{used_idxs(idx)}, 'R_Profile', 'RVIR', ...
        all_profiles{used_idxs(idx)}, 'T_Profile', 'TVIR');

%     plot_rg_params{1+2*(idx-1)} = all_profiles{used_idxs(idx)}.R_Profile'/all_profiles{used_idxs(idx)}.RVIR;
%     plot_rg_params{2+2*(idx-1)} = all_profiles{used_idxs(idx)}.RHOG_Profile'/all_profiles{used_idxs(idx)}.RHOVIR;
    [plot_rg_params{1+2*(idx-1)} plot_rg_params{2+2*(idx-1)}] = merge_multi_res_from_hash(...
        all_profiles{used_idxs(idx)}, 'R_Profile', 'RVIR', ...
        all_profiles{used_idxs(idx)}, 'RHOG_Profile', 'RHOVIR');
% 
%     plot_rt_params{1+2*(idx-1)} = all_profiles{used_idxs(idx)}.R_Profile'/all_profiles{used_idxs(idx)}.RVIR;
%     plot_rt_params{2+2*(idx-1)} = all_profiles{used_idxs(idx)}.RHOTOT_Profile'/all_profiles{used_idxs(idx)}.RHOVIR;
    [plot_rt_params{1+2*(idx-1)} plot_rt_params{2+2*(idx-1)}] = merge_multi_res_from_hash(...
        all_profiles{used_idxs(idx)}, 'R_Profile', 'RVIR', ...
        all_profiles{used_idxs(idx)}, 'RHOTOT_Profile', 'RHOVIR');
% 
%     plot_s_params{1+2*(idx-1)} = all_profiles{used_idxs(idx)}.R_Profile'/all_profiles{used_idxs(idx)}.RVIR;
%     plot_s_params{2+2*(idx-1)} = all_profiles{used_idxs(idx)}.Normalized_S_Profile';
    [plot_s_params{1+2*(idx-1)} plot_s_params{2+2*(idx-1)}] = merge_multi_res_from_hash(...
        all_profiles{used_idxs(idx)}, 'R_Profile', 'RVIR', ...
        all_profiles{used_idxs(idx)}, 'Normalized_S_Profile', '');

%     FVIR = G * all_profiles{used_idxs(idx)}.MVIR * Msun / ((all_profiles{used_idxs(idx)}.RVIR*MPc)^2)
    plot_turb_params{1+2*(idx-1)} = all_turbs{used_idxs(idx)}.bins_x_all;
    plot_turb_params{2+2*(idx-1)} = all_turbs{used_idxs(idx)}.bins;

%     %forces
%     plot_hse_params{1+2*(idx-1)} = all_profiles{used_idxs(idx)}.R2_Profile/all_profiles{used_idxs(idx)}.RVIR;
%     plot_hse_params{2+2*(idx-1)} = abs(all_profiles{used_idxs(idx)}.RIGHT./all_profiles{used_idxs(idx)}.LEFT);
    [r2_profile right_profile] = merge_multi_res_from_hash(...
        all_profiles{used_idxs(idx)}, 'R2_Profile', 'RVIR', ...
        all_profiles{used_idxs(idx)}, 'RIGHT', '');
    [r2_profile left_profile] = merge_multi_res_from_hash(...
        all_profiles{used_idxs(idx)}, 'R2_Profile', 'RVIR', ...
        all_profiles{used_idxs(idx)}, 'LEFT', '');
    plot_hse_params{1+2*(idx-1)} = r2_profile;
    plot_hse_params{2+2*(idx-1)} = right_profile ./ left_profile;
%     
%     %fluxes
% %     figure; plot(R_Profile/RVIR, [IN_FLOW_TEMP; OUT_FLOW_TEMP]', '.-'); title('Temperature profile of inflowing and outflowing fluxes');xlabel('r / R_{v}'); ylabel('Mean Temperature (K)'); legend('In Flow', 'Out Flow'); saveas(gcf, getresultsdir(sprintf('%s - InOut Flux Temperature(%d MPc).png',clustername, MPSec)));
% 
%     plot_flux_params{1+2*(idx-1)} = all_profiles{used_idxs(idx)}.R_Profile'/all_profiles{used_idxs(idx)}.RVIR;
%     plot_flux_params{2+2*(idx-1)} = all_fluxes{used_idxs(idx)}.M_DOT/abs(all_fluxes{used_idxs(idx)}.used_mdot);
    [plot_flux_params{1+2*(idx-1)} plot_flux_params{2+2*(idx-1)}] = merge_multi_res_from_hash(...
        all_profiles{used_idxs(idx)}, 'R_Profile', 'RVIR', ...
        all_fluxes{used_idxs(idx)}, 'M_DOT', 'used_mdot');
    plot_flux_params{2+2*(idx-1)} = plot_flux_params{2+2*(idx-1)} * (-1); % this compensates for the abs() above

end

%%
%     figure; plot(R_Profile/RVIR, [M_DOT/abs(used_mdot); Min_DOT/abs(used_mdot); Mout_DOT/abs(used_mdot)]', '.-'); 
% 
% line([min(R_Profile/RVIR) max(R_Profile/RVIR)], [0 0]); legend('Total', 'In', 'Out'); saveas(gcf, getresultsdir(sprintf('%s - Mdot InOut (%d MPc).png',clustername, MPSec)));
figure;
plot(plot_flux_params{:});
ylim([-5 2]);
grid on;
xlabel('r / R_{v}'); ylabel('Mdot / Mdot_{v}');
legend_h = legend(clusternames{used_idxs}, 'Location', 'NorthEast');set(legend_h, 'Interpreter', 'none');legend(legend_h);
title('Fluxes calculated in the 8Mpc cubes');
saveas(gcf,getresultsdir('all - fluxes.png'));


%%
figure;
loglog(plot_turb_params{:});
grid on;
xlabel('k (1/Mpc)');ylabel('E(k) (arbitrary)');
legend_h = legend(clusternames{used_idxs}, 'Location', 'SouthWest');set(legend_h, 'Interpreter', 'none');legend(legend_h);
title('Calculating the power spectrum of V^{2} in the 1Mpc cube');
set(gca,'YTick',10.^(10:25));
saveas(gcf,getresultsdir('all - Power spectrums.png'));

%%
figure;
plot(plot_hse_params{:});
grid on;
title('Hydrostatic Equilibrium for various clusters');
legend_h = legend(clusternames{used_idxs}, 'Location', 'NorthEast');set(legend_h, 'Interpreter', 'none');legend(legend_h);
ylabel('Pressure Term / Gravitational Term'); xlabel('r / R_{v}');
ylim([0 2]);xlim([0 2]);
saveas(gcf, getresultsdir('all - Hydrostatic Ratio (8Mpc).png'));

%%
figure;
loglog(plot_t_params{:}); 
legend_h = legend(clusternames{used_idxs}, 'Location', 'SouthWest');set(legend_h, 'Interpreter', 'none');legend(legend_h);
title('T profile for various clusters at z=0');
xlabel('r / R_{v}');ylabel('T / T_{v}');grid on;
saveas(gcf, getresultsdir('all - T profiles - 8Mpc.png'));

figure;
loglog(plot_rg_params{:}); 
legend_h = legend(clusternames{used_idxs}, 'Location', 'SouthWest');set(legend_h, 'Interpreter', 'none');legend(legend_h);
title('\rho_g profile for various clusters at z=0');
xlabel('r / R_{v}');ylabel('\rho_{g} / \rho_{v}');grid on;
saveas(gcf, getresultsdir('all - RHOG profiles - 8Mpc.png'));

figure;
loglog(plot_rt_params{:}); 
legend_h = legend(clusternames{used_idxs}, 'Location', 'SouthWest');set(legend_h, 'Interpreter', 'none');legend(legend_h);
title('\rho_{tot} profile for various clusters at z=0');
xlabel('r / R_{v}');ylabel('\rho / \rho_{v}');grid on;
saveas(gcf, getresultsdir('all - RHOTOT profiles - 8Mpc.png'));


figure;
loglog(plot_s_params{:}); 
legend_h = legend(clusternames{used_idxs}, 'Location', 'SouthEast');set(legend_h, 'Interpreter', 'none');legend(legend_h);
title('S profile for various clusters at z=0');
xlabel('r / R_{v}');ylabel('Gas Entropy K= (T/T_{v})/(\rho/\rho_{v})^{2/3})');grid on;
saveas(gcf, getresultsdir('all - S profiles - 8Mpc.png'));
