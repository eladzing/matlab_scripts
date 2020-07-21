%% Init
global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';
 

%% Load
Vrr = Vr(8);
RG = RHOG(8);
f2 = Vrr.*RG;
T2 = T(8);

%% Base Figures
figure(1); imshow(squeeze(sum(f2(end:-1:1,127:130,:),2)),[-2e14 2e14]); colormap(jet); colorbar();
figure(2); imshow(squeeze(sum(T2(end:-1:1,127:130,:),2)),[]); colormap(jet); colorbar();
draw_circle(1, 30);
draw_circle(1, 70);
draw_circle(2, 30);
draw_circle(2, 70);

%% FT Plots
[ff TT SS] = flux(2,30*4);
plot_FT(ff(:),TT(:),30*4,2);
figure; imshow(ff, 1e14*[-3 3]); colormap(jet); colorbar(); xlabel('Theta'); ylabel('Cos(Phi)'); title('Radial Flux at (R=~2MPSec)')

[ff TT SS] = flux(4,60*2);
plot_FT(ff(:),TT(:),60*2,4);
figure; imshow(ff, 1e14*[-2 2]); colormap(jet); colorbar(); xlabel('Theta'); ylabel('Cos(Phi)'); title('Radial Flux at (R=~1MPSec)')

% [ff TT SS] = flux(8,70);
% plot_FT(ff(:),TT(:),70,8);
% figure; imshow(ff, 1e14*[-2 2]); colormap(jet);

%%
% for R = [120:-10:10]
% [ff TT SS] = flux(8,R); figure(4); 
% [nff ticks] = normalize_colormap(ff);
% imshow(nff); colormap(jet); % colorbar();
% pause(1/100);
% end

%%
IN_FLOW_TEMP = [];
OUT_FLOW_TEMP = [];
Min_DOT  = [];
Mout_DOT = [];
Mcold_DOT = [];
M_DOT    = [];
RR = [];
MPSec = 2;
for R = [120:-5:5]
    R
    RR = [RR R*(MPSec*1000)/256];
    [ff TT SS] = flux(MPSec,R); 
    IN_FLOW_TEMP  = [IN_FLOW_TEMP  sum(TT(ff<0).*ff(ff<0))/sum(ff(ff<0))];
    OUT_FLOW_TEMP = [OUT_FLOW_TEMP sum(TT(ff>0).*ff(ff>0))/sum(ff(ff>0))];
    Min_DOT       = [Min_DOT  sum(ff(ff<0))];
    Mout_DOT      = [Mout_DOT sum(ff(ff>0))];
    M_DOT         = [M_DOT    sum(ff(:))];
    
    cold_mask = bitand(TT<1.5e7,SS<0.135);
    Mcold_DOT     = [Mcold_DOT   sum(ff(cold_mask))];
end

figure; plot(RR/1000, [IN_FLOW_TEMP; OUT_FLOW_TEMP]', '*-'); xlabel('R (MPSec)'); ylabel('Mean Temperature (K)'); legend('In Flow', 'Out Flow')
figure; plot(RR/1000, abs([Min_DOT; Mout_DOT])', '*-'); xlabel('R (MPSec)'); ylabel('M dot'); legend('In', 'Out')
figure; plot(RR/1000, [M_DOT]', '*-'); xlabel('R (MPSec)'); ylabel('M dot')

%% "AITOFF"
