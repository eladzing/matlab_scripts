%%
global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

RR8 = RHOG(8);
TT8 = T(8);
SS8 = S(8);
VR8 = load('../data/Vr_sphere_8','res'); VR8 = VR8.res;
T8 = load('../data/T_sphere_8','res'); T8 = T8.res;
S8 = load('../data/S_sphere_8','res'); S8 = S8.res;
ds8 = ds_sphere(8);
RG8 = load('../data/RHOG_sphere_8','res'); RG8 = RG8.res;

%%
figure(1); [range vals] = hist_density((T8(64:end-74,:,:)),RG8(64:end-74,:,:).*ds8(64:end-74,:,:),100,0,0.6e7); bar(range,vals); saveas(gcf, getresultsdir('BiModality - T non-logarithmic.png'));
figure(2); imshow(squeeze(sum(TT8(end-74:-1:64,64:end-74,64:end-74)<5.7e5,2)), [], 'Border', 'tight'); colormap(jet); colorbar; saveas(gcf, getresultsdir('BiModality - T less than 5.7e5.png'));

figure(3); [range vals] = hist_density((log10(T8(64:end-74,:,:))),RG8(64:end-74,:,:).*ds8(64:end-74,:,:),150); bar(range,log10(vals)); saveas(gcf, getresultsdir('BiModality - T Logarithmic.png'));
figure(4); imshow(squeeze(sum(TT8(end-74:-1:64,64:end-74,64:end-74)<10^(3.83),2)), [], 'Border', 'tight'); colormap(jet); colorbar; saveas(gcf, getresultsdir('BiModality - T less than 1e3.83 Logarithmic.png'));

%%
figure(1); [range vals] = hist_density((S8(64:end-74,:,:)),RG8(64:end-74,:,:).*ds8(64:end-74,:,:),100,0,0.35); bar(range,vals); saveas(gcf, getresultsdir('BiModality - S non-logarithmic.png'));
figure(2); [range vals] = hist_density(log10(S8(64:end-74,:,:)),RG8(64:end-74,:,:).*ds8(64:end-74,:,:),100); bar(range,log10(vals)); saveas(gcf, getresultsdir('BiModality - S logarithmic.png'));

figure;
subplot(1,3,1);imshow(squeeze(sum(SS8(end:-1:1,127:131,:)<0.135,2)), [], 'Border', 'tight'); colormap(jet);
subplot(1,3,2);imshow(squeeze(sum(SS8(end:-1:1,:,127:131)<0.135,3)), [], 'Border', 'tight'); colormap(jet);
subplot(1,3,3);imshow(squeeze(sum(SS8(127:131,:,:)<0.135,1)), [], 'Border', 'tight'); colormap(jet);

%%
figure(1); [range vals] = hist_density((RG8(64:end-74,:,:)),RG8(64:end-74,:,:).*ds8(64:end-74,:,:),100); bar(range,vals); saveas(gcf, getresultsdir('BiModality - RhoG non-logarithmic.png'));
figure(1); [range vals] = hist_density(log(RG8(64:end-74,:,:)),RG8(64:end-74,:,:).*ds8(64:end-74,:,:),100); bar(range,vals); saveas(gcf, getresultsdir('BiModality - RhoG logarithmic.png'));
figure(1); imshow(squeeze(sum(RR8(end:-1:1,:,:)>1e12,2)), [], 'Border', 'tight'); colormap(jet);saveas(gcf, getresultsdir('BiModality - RhoG more than 1e12 logarithmic.png'));

%%
figure(1); [range vals] = hist_density((VR8(64:end-74,:,:)),RG8(64:end-74,:,:).*ds8(64:end-74,:,:),100); bar(range,vals); saveas(gcf, getresultsdir('BiModality - Vr non-logarithmic.png'));




%%%%%%%%%%%%%%%%%
mask1 = bitand(SS8(end:-1:1,127:131,:)<0.135, TT8(end:-1:1,127:131,:)<1.5e7);
mask2 = bitand(SS8(end:-1:1,:,127:131)<0.135, TT8(end:-1:1,:,127:131)<1.5e7);
mask3 = bitand(SS8(127:131,:,:)<0.135, TT8(127:131,:,:)<1.5e7);
figure;
subplot(1,3,1);imshow(squeeze(sum(mask1,2)), [], 'Border', 'tight'); colormap(jet);
subplot(1,3,2);imshow(squeeze(sum(mask2,3)), [], 'Border', 'tight'); colormap(jet);
subplot(1,3,3);imshow(squeeze(sum(mask3,1)), [], 'Border', 'tight'); colormap(jet);
saveas(gcf, getresultsdir('BiModality - S less than 0.135 AND T less than 1.5e7.png'));