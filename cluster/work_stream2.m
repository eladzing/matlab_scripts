global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

MPSec = 8;
%%%b = S(MPSec);
%%%mr = squeeze(sum(b(:,:,127:130),3));

%%
[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);
[Vrr cosa] = Vr_mesh(MPSec);
Vabs = sqrt(Vxx.^2 + Vyy.^2 + Vzz.^2);

%%
%mask = bitand(cosa<0, Vabs>200);
mask = (cosa<0);
Vxx(~mask) = 0;
Vyy(~mask) = 0;
Vzz(~mask) = 0;


% %% 3
% Vxx_slice = squeeze(sum(Vxx(:,:,127:130),3));
% Vyy_slice = squeeze(sum(Vyy(:,:,127:130),3));
% Vzz_slice = squeeze(sum(Vzz(:,:,127:130),3));
% 
% DILUTE = 4;
% figure; quiver(Vyy_slice(1:DILUTE:end,1:DILUTE:end),Vxx_slice(1:DILUTE:end,1:DILUTE:end),2.5)
% 
% %% 2
% Vxx_slice = squeeze(sum(Vxx(:,127:130,:),2));
% Vyy_slice = squeeze(sum(Vyy(:,127:130,:),2));
% Vzz_slice = squeeze(sum(Vzz(:,127:130,:),2));
% 
% DILUTE = 4;
% figure; quiver(Vzz_slice(1:DILUTE:end,1:DILUTE:end),Vxx_slice(1:DILUTE:end,1:DILUTE:end),2.5)
% 
% %% 1
% Vxx_slice = squeeze(sum(Vxx(127:130,:,:),1));
% Vyy_slice = squeeze(sum(Vyy(127:130,:,:),1));
% Vzz_slice = squeeze(sum(Vzz(127:130,:,:),1));
% 
% DILUTE = 4;
% figure; quiver(Vzz_slice(1:DILUTE:end,1:DILUTE:end),Vyy_slice(1:DILUTE:end,1:DILUTE:end),2.5)
% 
% 
% %%
% DILUTE = 1;
% [startX1 startY1] = meshgrid([1 256], [1:DILUTE:256]);
% [startX2 startY2] = meshgrid([1:DILUTE:256], [1 256]);
% vecs = stream2(Vyy_slice, Vxx_slice, [startX1(:) startX2(:)], [startY1(:) startY2(:)]);
% figure; streamline(vecs);

%%
DILUTE = 5;
% [startX1 startY1 startZ1] = meshgrid([1 256], [1:DILUTE:256], [1:DILUTE:256]);
% [startX2 startY2 startZ2] = meshgrid([1:DILUTE:256], [1 256], [1:DILUTE:256]);
% [startX3 startY3 startZ3] = meshgrid([1:DILUTE:256], [1:DILUTE:256], [1 256]);

[startY1 startX1 startZ1] = meshgrid([1:DILUTE:256], [1 256], [1:DILUTE:256]);
[startY2 startX2 startZ2] = meshgrid([1 256], [1:DILUTE:256], [1:DILUTE:256]);
[startY3 startX3 startZ3] = meshgrid([1:DILUTE:256], [1:DILUTE:256], [1 256]);


%%
[meshY, meshX, meshZ] = meshgrid(1:size(Vyy,1), 1:size(Vxx,2), 1:size(Vzz,3));
%verts = stream3(meshX, meshY, meshZ, Vxx, Vyy, Vzz, [startX1(:) startX2(:) startX3(:)], [startY1(:) startY2(:) startY3(:)], [startZ1(:) startZ2(:) startZ3(:)]);
verts = stream3(meshY, meshX, meshZ, Vyy, Vxx, Vzz, [startY1(:) startY2(:) startY3(:)], [startX1(:) startX2(:) startX3(:)], [startZ1(:) startZ2(:) startZ3(:)]);
mask = f(verts,40);
mask = logical(mask);
