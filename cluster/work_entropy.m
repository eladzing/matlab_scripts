global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

MPSec = 8;
b = S(MPSec);
%%%slice = transpose(squeeze(sum(b(127:130,:,end:-1:1),1)));
%%%imshow(log(slice), []); colormap(jet);
%%%pause;

% TODO % prepend black->blue and append red->white

%%% TRY1: based on jet colormap (UGLY!)
% post = ([0:127]'./127)*([0 1 1]);
% post(:,1) = 1;
% pre = ([0:127]'./127)*([0 0 1]);
% 
% jet2 = zeros(512,3);
% jet2(1:128,:)         = pre;
% jet2(129:129+256-1,:) = jet(256);
% jet2(129+256:end,:)   = post;
% 
% imshow(log(slice), [-1.5 4.5]); colormap(jet2);



%%%1111 fix%%%
Vx = read_cube(sprintf(FILE_FORMAT, 'vx', MPSec));
Vy = read_cube(sprintf(FILE_FORMAT, 'vy', MPSec));
Vz = read_cube(sprintf(FILE_FORMAT, 'vz', MPSec));
figure(1); slice = transpose(squeeze(sum(b(:,end:-1:1,127:130),3))); imshow(log(slice(end:-1:1,:)), []); colormap(jet);
Vx = transpose(squeeze(sum(Vx.data(:,end:-1:1,127:130),3)));
Vy = transpose(squeeze(sum(Vy.data(:,end:-1:1,127:130),3)));
Vz = transpose(squeeze(sum(Vz.data(:,end:-1:1,127:130),3)));
figure(2); quiver(Vx(1:8:end,1:8:end),(-1)*Vy(1:8:end,1:8:end),1.5)
%%%%%%

%%%2222%%%
Vx = read_cube(sprintf(FILE_FORMAT, 'vx', MPSec));
Vy = read_cube(sprintf(FILE_FORMAT, 'vy', MPSec));
Vz = read_cube(sprintf(FILE_FORMAT, 'vz', MPSec));
figure(1); slice = transpose(squeeze(sum(b(127:130,:,end:-1:1),1))); imshow(log(slice), []); colormap(jet);
Vx = transpose(squeeze(sum(Vx.data(127:130,:,end:-1:1),1)));
Vy = transpose(squeeze(sum(Vy.data(127:130,:,end:-1:1),1)));
Vz = transpose(squeeze(sum(Vz.data(127:130,:,end:-1:1),1)));
figure(2); quiver(Vy(end:-8:1,1:8:end),Vz(end:-8:1,1:8:end),3)
%%%%%%%

%%%3333%%%%
Vx = read_cube(sprintf(FILE_FORMAT, 'vx', MPSec));
Vy = read_cube(sprintf(FILE_FORMAT, 'vy', MPSec));
Vz = read_cube(sprintf(FILE_FORMAT, 'vz', MPSec));
figure(1); slice = transpose(squeeze(sum(b(:,127:130,end:-1:1),2))); imshow(log(slice), []); colormap(jet);
Vx = transpose(squeeze(sum(Vx.data(:,127:130,end:-1:1),2)));
Vy = transpose(squeeze(sum(Vy.data(:,127:130,end:-1:1),2)));
Vz = transpose(squeeze(sum(Vz.data(:,127:130,end:-1:1),2)));
figure(2); quiver(Vx(end:-8:1,1:8:end),Vz(end:-8:1,1:8:end),2)
%%%%%%%


%%%%%3333 fix%%%%%%%%%%%%%%%%%%%%55
Vx = read_cube(sprintf(FILE_FORMAT, 'vx', MPSec));
Vy = read_cube(sprintf(FILE_FORMAT, 'vy', MPSec));
Vz = read_cube(sprintf(FILE_FORMAT, 'vz', MPSec));
figure(1); slice = transpose(squeeze(sum(b(:,127:130,:),2))); imshow(log(slice(end:-1:1,:)), []); colormap(jet);
Vx = transpose(squeeze(sum(Vx.data(:,127:130,:),2)));
Vy = transpose(squeeze(sum(Vy.data(:,127:130,:),2)));
Vz = transpose(squeeze(sum(Vz.data(:,127:130,:),2)));
figure(2); quiver(Vx(1:8:end,1:8:end),Vz(1:8:end,1:8:end),2)



%%%2222 fix%%%
Vx = read_cube(sprintf(FILE_FORMAT, 'vx', MPSec));
Vy = read_cube(sprintf(FILE_FORMAT, 'vy', MPSec));
Vz = read_cube(sprintf(FILE_FORMAT, 'vz', MPSec));
figure(1); slice = transpose(squeeze(sum(b(127:130,end:-1:1,:),1))); imshow(log(slice(end:-1:1,:)), []); colormap(jet);
Vx = transpose(squeeze(sum(Vx.data(127:130,end:-1:1,:),1)));
Vy = transpose(squeeze(sum(Vy.data(127:130,end:-1:1,:),1)));
Vz = transpose(squeeze(sum(Vz.data(127:130,end:-1:1,:),1)));
figure(2); quiver(Vy(1:8:end,1:8:end),Vz(1:8:end,1:8:end),3)
%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% 3333 fix (no transpose)
figure(1); slice = squeeze(sum(b(:,127:130,:),2)); imshow(log(slice(end:-1:1,:)), []); colormap(jet);
Vxx = Vx(MPSec); Vxx = squeeze(sum(Vxx(:,127:130,:),2));
Vyy = Vx(MPSec); Vyy = squeeze(sum(Vyy(:,127:130,:),2));
Vzz = Vz(MPSec); Vzz = squeeze(sum(Vzz(:,127:130,:),2));
figure(2); quiver(Vzz(1:8:end,1:8:end),Vxx(1:8:end,1:8:end),2)

%% 1111 fix (no transpose)
figure(1); slice = squeeze(sum(b(:,:,127:130),3)); imshow(log(slice(end:-1:1,:)), []); colormap(jet);
Vxx = Vx(MPSec); Vxx = squeeze(sum(Vxx(:,:,127:130),3));
Vyy = Vy(MPSec); Vyy = squeeze(sum(Vyy(:,:,127:130),3));
Vzz = Vz(MPSec); Vzz = squeeze(sum(Vzz(:,:,127:130),3));
figure(2); quiver(Vyy(1:8:end,1:8:end),Vxx(1:8:end,1:8:end),1.5)

%% 2222 fix (no transpose)
figure(1); slice = squeeze(sum(b(127:130,:,:),1)); imshow(log(slice(end:-1:1,:)), []); colormap(jet);
Vxx = Vx(MPSec); Vxx = squeeze(sum(Vxx(127:130,:,:),1));
Vyy = Vy(MPSec); Vyy = squeeze(sum(Vyy(127:130,:,:),1));
Vzz = Vz(MPSec); Vzz = squeeze(sum(Vzz(127:130,:,:),1));
figure(2); quiver(Vzz(1:8:end,1:8:end),Vyy(1:8:end,1:8:end),3)