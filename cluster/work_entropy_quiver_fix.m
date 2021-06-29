global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

%for MPSec = [1 2 4 8]
MPSec = 8;
%mr = multi_res('S(%d)', 'sum(%s(:,127:130,:),2)');
b = S(MPSec);
mr = squeeze(sum(b(:,127:130,:),2));
%mr = imresize(mr, 8);
%Vx = read_cube(sprintf(FILE_FORMAT, 'vx', MPSec));
%Vy = read_cube(sprintf(FILE_FORMAT, 'vy', MPSec));


%%%%%%%%%%%%%%%%%5
%figure(1); slice = squeeze(sum(b(:,127:130,:),2)); imshow(log(slice(end:-1:1,:)), []); colormap(jet);
% Vxx = Vx(MPSec); Vxx = squeeze(sum(Vxx(:,127:130,:),2));
% Vyy = Vx(MPSec); Vyy = squeeze(sum(Vyy(:,127:130,:),2));
% Vzz = Vz(MPSec); Vzz = squeeze(sum(Vzz(:,127:130,:),2));

[Vxx Vyy Vzz VcmX VcmY VcmZ] = V_Vcm(MPSec);
Vxx = squeeze(sum(Vxx(:,127:130,:),2));
Vyy = squeeze(sum(Vyy(:,127:130,:),2));
Vzz = squeeze(sum(Vzz(:,127:130,:),2));

%figure(2); quiver(Vzz(1:8:end,1:8:end),Vxx(1:8:end,1:8:end),2)
%%%%%%%%%%%%%%%%%5



%%

%%%%%%%% GRAD TO NO GRAD
%[FX FY] = gradient(mr);
%b = sqrt(FX.^2 + FY.^2);
%b = log(1+b(end:-1:1,:));

b = mr;%%%%%%%% GRAD TO NO GRAD

%%imshow(log(b(end:-1:1,:)), []); colormap(jet); pause;
%b(b<0) = 0;

% Strech
b = log(b(end:-1:1,:));   %%%%%%%% GRAD TO NO GRAD
b = imresize(b, 8);
b = b - min(b(:));
b = b ./ min(max(b(:)));  %%%%%%%%%%%%0.1*MPSec, 
b = ind2rgb(int32(1+b*(2^8-1)), jet); %2.^16, hot
%imshow(b); pause;

imshow(b);
%pause

DILUTE = 4;
%quiver(Vx.data(1:8:end,1:8:end,128),Vy.data(1:8:end,1:8:end,128),1.5)
quiver(Vzz(1:DILUTE:end,1:DILUTE:end),Vxx(1:DILUTE:end,1:DILUTE:end),2.5)
xlim([3 256/DILUTE]); ylim([2 256/DILUTE-1]);
f = getframe();
f = f.cdata;

% white -> black, blue -> white
f = rgb2gray(f);
f(f==255) = 0;
f(f==29) = 255;
[X M] = gray2ind(f);
f = im2double(ind2rgb(X,M));

% resize to fit
f = imresize(f, max(size(b))/min(size(f,1), size(f,2)));
f = f(1:size(b,1), 1:size(b,2), :);

whos f
whos b

b = b + f;
b(b>1) = 1;
imshow(b);
imwrite(b, sprintf('entropy_quiver_fix_%d.jpg',MPSec));

