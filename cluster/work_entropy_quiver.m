global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

MPSec = 8;
mr = multi_res('S(%d)', '(:,:,128)');
Vx = read_cube(sprintf(FILE_FORMAT, 'vx', MPSec));
Vy = read_cube(sprintf(FILE_FORMAT, 'vy', MPSec));


%%
b = mr;
b(b<0) = 0;
b = log10(1+b);
b = b - min(b(:));
b = b ./ max(b(:));
b = ind2rgb(int32(1+b*(2^8-1)), hot); %2.^16
imshow(b);
pause

quiver(Vx.data(1:8:end,1:8:end,128),Vy.data(1:8:end,1:8:end,128),1.5)
xlim([0 33]); ylim([0 33]);
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
imwrite(b, 'entropy_quiver.jpg');