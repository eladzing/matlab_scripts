function batch_turbulence(halopath, clustername)
%%
MPSec = 1;
[vx, vy, vz, v1,v2,v3] = V_Vcm(MPSec);
vel_abs = (vx.^2+vy.^2+vz.^2);

%%
%vel_abs = vel_abs.*RHOG(1);
%%
[bins_x_all   bins] = calc_kolmogorov(vel_abs);
[bins_x_small bins_small] = calc_kolmogorov(vel_abs(65:192,65:192,65:192));

[bins_x bins1] = calc_kolmogorov(vel_abs(1:128,   1:128,   1:128));

[bins_x bins2] = calc_kolmogorov(vel_abs(129:end, 1:128,   1:128));
[bins_x bins3] = calc_kolmogorov(vel_abs(1:128,   129:end, 1:128));
[bins_x bins4] = calc_kolmogorov(vel_abs(1:128,   1:128,   129:end));

[bins_x bins5] = calc_kolmogorov(vel_abs(129:end, 129:end, 1:128));
[bins_x bins6] = calc_kolmogorov(vel_abs(129:end, 1:128,   129:end));
[bins_x bins7] = calc_kolmogorov(vel_abs(1:128,   129:end, 129:end));

[bins_x bins8] = calc_kolmogorov(vel_abs(129:end,129:end,129:end));



%%
% loglog(bins_x_all, bins, 'o-', 10:80, (10:80).^(-2.5).*(10^(12.7)),'--', bins_x, [bins1;bins2;bins3;bins4;bins5;bins6;bins7;bins8]);
% grid on
% %loglog(bins_x_all, bins, 'o-', 10:80, (10:80).^(-2.5).*(10^(28)),'--', bins_x, [bins1;bins2;bins3;bins4;bins5;bins6;bins7;bins8]);;
% xlabel('k (1/Mpc)');ylabel('E(k) (arbitrary)');
% legend('entire 1Mpc cube','reference slope of -2.5','subcube1','subcube2','subcube3','subcube4','subcube5','subcube6','subcube7','subcube8','Location', 'SouthWest')
% title('Calculating the power spectrum of V^{2}')
% saveas(gcf,getresultsdir('Calculating the power spectrum of V2.png'))

figure;
loglog(bins_x_all, bins, 'o-', bins_x_small, bins_small, 'o-', bins_x, [bins1;bins2;bins3;bins4;bins5;bins6;bins7;bins8], '--');
grid on
xlabel('k (1/Mpc)');ylabel('E(k) (arbitrary)');
legend('entire 1Mpc cube','middle 0.5 Mpc cube','subcube1','subcube2','subcube3','subcube4','subcube5','subcube6','subcube7','subcube8','Location', 'SouthWest')
title('Calculating the power spectrum of V^{2}')
saveas(gcf,getresultsdir(sprintf('%s - Calculating the power spectrum of V2.png', clustername)));

save(sprintf('%s/turbulence',halopath), 'bins_x_all', 'bins', 'bins_x_small', 'bins_small', 'bins_x', 'bins1', 'bins2', 'bins3', 'bins4', 'bins5', 'bins6', 'bins7', 'bins8');