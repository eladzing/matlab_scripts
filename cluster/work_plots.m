iptsetpref('ImshowBorder','tight')
global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';


%%
figure();
TT = T(8); subplot(2,2,1); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[6 8]); colormap(jet); colorbar(); title('8 MPSec');
TT = T(4); subplot(2,2,2); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[6 8]); colormap(jet); colorbar(); title('4 MPSec');
TT = T(2); subplot(2,2,3); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[6 8]); colormap(jet); colorbar(); title('2 MPSec');
TT = T(1); subplot(2,2,4); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[6 8]); colormap(jet); colorbar(); title('1 MPSec');

%%
figure();
%title('log(\rho_{gas}) Multiresolution (8,4,2,1 Mpc/h)')
TT = RHOG(8); subplot(2,2,1); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[9 15]); colormap(jet); colorbar(); title('8 MPc/h');
TT = RHOG(4); subplot(2,2,2); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[9 15]); colormap(jet); colorbar(); title('4 MPc/h');
TT = RHOG(2); subplot(2,2,3); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[9 15]); colormap(jet); colorbar(); title('2 MPc/h');
TT = RHOG(1); subplot(2,2,4); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[9 15]); colormap(jet); colorbar(); title('1 MPc/h');
saveas(gcf, getresultsdir('Rhog Profile (multi resolution).png'))

%%
figure();
%title('log(\rho_{gas}) Multiresolution (8,4,2,1 Mpc/h)')
TT = RHOTOT(8); subplot(2,2,1); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[8 16]); colormap(jet); colorbar(); title('8 MPc/h');
TT = RHOTOT(4); subplot(2,2,2); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[8 16]); colormap(jet); colorbar(); title('4 MPc/h');
TT = RHOTOT(2); subplot(2,2,3); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[8 16]); colormap(jet); colorbar(); title('2 MPc/h');
TT = RHOTOT(1); subplot(2,2,4); imshow(squeeze(sum(log10(TT(end:-1:1,127:130,:)),2)/4),[8 16]); colormap(jet); colorbar(); title('1 MPc/h');
saveas(gcf, getresultsdir('Rhotot Profile (multi resolution).png'))

%%
figure();
TT = Vr(8); subplot(2,2,1); imshow(squeeze(sum((TT(end:-1:1,127:130,:)),2)/4),[-1500 750]); colormap(jet); colorbar(); title('8 MPc/h');
TT = Vr(4); subplot(2,2,2); imshow(squeeze(sum((TT(end:-1:1,127:130,:)),2)/4),[-1500 750]); colormap(jet); colorbar(); title('4 MPc/h');
TT = Vr(2); subplot(2,2,3); imshow(squeeze(sum((TT(end:-1:1,127:130,:)),2)/4),[-1500 750]); colormap(jet); colorbar(); title('2 MPc/h');
TT = Vr(1); subplot(2,2,4); imshow(squeeze(sum((TT(end:-1:1,127:130,:)),2)/4),[-1500 750]); colormap(jet); colorbar(); title('1 MPc/h');
saveas(gcf, getresultsdir('V_radial Profile (multi resolution).png'))