global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

MPSec = 1;
Vx = read_cube(sprintf(FILE_FORMAT, 'vx', MPSec));
Vy = read_cube(sprintf(FILE_FORMAT, 'vy', MPSec));
Vz = read_cube(sprintf(FILE_FORMAT, 'vz', MPSec));
V = sqrt((Vx.data).^2 + (Vy.data).^2 + (Vz.data).^2);
%clear Vx Vy Vz

Vf = fftshift(fftn(V));
imshow(log(1+abs(Vf(:,:,128))),[]); colormap(hot);

%Vf = fftshift(fft2(V(:,:,128)));

%imshow(V(:,:,128),[]);colormap(hot);
%pause;
%imshow(log(1+abs(Vf)),[]); colormap(hot);

Vf2 = abs(Vf.^2);
vals = [];

for margin=1:128
    vals = [vals sum(sum(sum(Vf2(margin:end-margin+1, margin:end-margin+1, margin:end-margin+1))))];
end

plot((128:-1:1).^(-2/3), vals, '*-')