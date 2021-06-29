global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

R = RHOG(8);

global xmesh;
global ymesh;
global zmesh;
[xmesh ymesh zmesh] = meshgrid(1:size(R,1),1:size(R,2),1:size(R,3));


MOVIE = avifile('anim.avi','fps',2);%,'COMPRESSION','Cinepak');

%[min(log(1+im(:))) max(log(1+im(:)))]

minmax = [26.7325   35.8837]
%figure;
for a = 0:pi/100:2*pi
    im = zeros(256,256);
    for depth = 1:1:256
        im = im + warp_mat3d_fast(R, angle2dcm4(0,a,0), 1, depth);
    end
    im(im==0)=NaN;
    
    imshow(log(1+im),minmax); colormap(hot);
    frm = getframe();
    MOVIE = addframe(MOVIE,frm.cdata);
    pause(1/1000);
end


MOVIE = close(MOVIE);
