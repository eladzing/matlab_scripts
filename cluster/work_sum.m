global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

R = RHOG(8);

global xmesh;
global ymesh;
global zmesh;
[xmesh ymesh zmesh] = meshgrid(1:size(R,1),1:size(R,2),1:size(R,3));


%[min(log(1+im(:))) max(log(1+im(:)))]

minmax = [26.7325   35.8837]
figure;
for a = pi/2:pi/2:pi/2%2*pi
    im = zeros(256,256);
    for depth = 1:10:256
        im = im + warp_mat3d_fast(R, angle2dcm4(0,a,0), 1, depth);
    end
    im(im==0)=NaN;
    
    imshow(log(1+im),minmax); colormap(hot);
    pause(1/100);
end