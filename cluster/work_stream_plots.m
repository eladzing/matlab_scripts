function work_stream_plots(verts)

load profiles8
T8 = T(8);
S8 = S(8);
RG8 = RHOG(8);

%%
figure;
[T_results T_stds T_min T_max] = stream_profile_bins(verts, 40, T8, length(verts));
XX = [1:length(T_results)]/100;

%% Variance
loglog(XX, T_results, XX, T_results+T_stds, 'r.', XX, T_results-T_stds, 'r.', R_Profile/2.2, T_Profile)
xlabel('r/R_{v}')
ylabel('T')
legend('Average temperature of stream lines', 'STD bars','STD bars', 'average temperature of the cube', 'Location', 'SouthWest');
% ylim([min([T_results,T_Profile]) max([T_results,T_Profile])])
title('Temperature of selected streamlines vs. cube')
saveas(gcf, getresultsdir('Teperature of streamlines vs. cube (with stds).png'))


% figure;
% %%
% loglog(XX, T_results, XX, T_min, '.', XX, T_max, '.', R_Profile/2.2, T_Profile)
% xlabel('r/R_{v}')
% ylabel('T')
% legend('Average temperature of stream lines', 'min temperature of streamlines','max temperature of streamlines', 'average temperature of the cube', 'Location', 'SouthWest');
% ylim([min([T_results,T_Profile]) max([T_results,T_Profile])])
% title('Temperature of streamlines vs. cube')
% saveas(gcf, getresultsdir('Teperature of streamlines vs. cube.png'))

%%
 [S_results S_stds S_min S_max] = stream_profile_bins(verts, 40, S8, length(verts));
% %%
figure;
XX = [1:length(S_results)]/100;
loglog(XX, S_results, XX, S_results+S_stds, 'r.', XX, S_results-S_stds, 'r.', R_Profile/2.2, S_Profile)
xlabel('r/R_{v}')
ylabel('S')
legend('Average entropy of stream lines', 'STD bars','STD bars', 'average entropy of the cube', 'Location', 'NorthWest');
ylim([min([S_results,S_Profile]) max([S_results,S_Profile])])
title('Entropy of selected streamlines vs. cube')
saveas(gcf, getresultsdir('Entropy of streamlines vs. cube (with stds).png'))


%%
[RHOG_results RHOG_stds RHOG_min RHOG_max] = stream_profile_bins(verts, 40, RG8, length(verts));
% %%
figure;
XX = [1:length(RHOG_results)]/100;
loglog(XX, RHOG_results, XX, RHOG_results+RHOG_stds, 'r.', XX, RHOG_results-RHOG_stds, 'r.', R_Profile/2.2, RHOG_Profile)
xlabel('r/R_{v}')
ylabel('\rho')
legend('Average density of stream lines', 'STD bars','STD bars', 'average density of the cube', 'Location', 'NorthWest');
ylim([min([RHOG_results,RHOG_Profile]) max([RHOG_results,RHOG_Profile])])
title('Density of selected streamlines vs. cube')
saveas(gcf, getresultsdir('Density of streamlines vs. cube (with stds).png'))
