%% Init
global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

%%
[ff TT SS] = flux(2,30*4);

%%
%mask = bitand(ff<1e15,ff>-1e15); TTT = TT(mask); fff = ff(mask);
TTT = TT(:); fff = ff(:);

%%
%plot(log(TTT), fff, '.')
figure;
scattercloud(log10(TTT), fff, 50, 1, 'k+', jet(2^16), 100)
title('F-T Scatter with density (R~1MPc)');
xlabel('log(T)');
ylabel('Flux (Arbitrary Units)');
%saveas(gcf, getresultsdir('F-T Scatter with Density at R=0.94MPc.png'));

%%
[ff TT SS] = flux(4,60*2);

%%
mask = bitand(ff<1e15,ff>-1e15);

TTT = TT(mask);
fff = ff(mask);

%plot(log(TTT), fff, '.')
figure;
scattercloud(log10(TTT), fff, 50, 1, 'k+', jet(2^16), 100)
line([min(log10(TTT)) max(log10(TTT))], [0 0]);
title('F-T Scatter with density (R~2MPc)');
xlabel('log(T)');
ylabel('Flux (Arbitrary Units)');
%saveas(gcf, getresultsdir('F-T Scatter with Density at R=1.88MPc.png'));