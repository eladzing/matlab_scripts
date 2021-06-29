%%
MPSec = 1;
[vx, vy, vz, v1,v2,v3] = V_Vcm(MPSec);
vel_abs = (vx.^2+vy.^2+vz.^2);

%%
%vel_abs = vel_abs.*RHOG(1);
%%
[bins_x_all bins] = calc_kolmogorov(vel_abs);

[bins_x bins1] = calc_kolmogorov(vel_abs(1:128,   1:128,   1:128));

[bins_x bins2] = calc_kolmogorov(vel_abs(129:end, 1:128,   1:128));
[bins_x bins3] = calc_kolmogorov(vel_abs(1:128,   129:end, 1:128));
[bins_x bins4] = calc_kolmogorov(vel_abs(1:128,   1:128,   129:end));

[bins_x bins5] = calc_kolmogorov(vel_abs(129:end, 129:end, 1:128));
[bins_x bins6] = calc_kolmogorov(vel_abs(129:end, 1:128,   129:end));
[bins_x bins7] = calc_kolmogorov(vel_abs(1:128,   129:end, 129:end));

[bins_x bins8] = calc_kolmogorov(vel_abs(129:end,129:end,129:end));

%%
loglog(bins_x_all, bins, 'o-', 10:80, (10:80).^(-2.5).*(10^(12.7)),'--', bins_x, [bins1;bins2;bins3;bins4;bins5;bins6;bins7;bins8]);;
%loglog(bins_x_all, bins, 'o-', 10:80, (10:80).^(-2.5).*(10^(28)),'--', bins_x, [bins1;bins2;bins3;bins4;bins5;bins6;bins7;bins8]);;
xlabel('k (1/Mpc)');ylabel('E(k) (arbitrary)');
legend('entire 1Mpc cube','reference slope of -2.5 (kolnmogorov)','subcube1','subcube2','subcube3','subcube4','subcube5','subcube6','subcube7','subcube8','Location', 'SouthWest')
title('Calculating the power spectrum of V^{2}')
saveas(gcf,getresultsdir('Calculating the power spectrum of V2.png'))
%loglog(bins_x, [bins1/bins1(1);bins2/bins2(1);bins3/bins3(1);bins4/bins4(1);bins5/bins5(1);bins6/bins6(1);bins7/bins7(1);bins8/bins8(1)]);legend('show')

%% Calculate energy of higher frequencies
%0.7 is the dk?
%energy = sum(bins(1:221).*4.*pi.*(bins_x_all(1:221).^2))*0.7;
%energy2 = mean(vel_abs(:))*((1/h)^3)

