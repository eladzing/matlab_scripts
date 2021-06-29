env
%%
work_stream2

%%
SS = S(8);
TT= T(8);
new_verts = filter_verts(verts, 40, SS, 0.1,0.5);
new_verts2 = filter_verts(new_verts, 40, log10(TT), 0.5,inf);

clear SS TT

mask = f(new_verts2,40);
%mask = logical(mask);


%%
global xmesh;
global ymesh;
global zmesh;
[xmesh ymesh zmesh] = meshgrid(1:256,1:256,1:256);

%%

% MOVIE = avifile('anim.avi','fps',2);%,'COMPRESSION','Cinepak');

%[min(log(1+im(:))) max(log(1+im(:)))]

minval = 0;
maxval = 16;

counter = 0;
for a = 0:pi/100:2*pi
    a
    im = zeros(256,256);
    for depth = 1:1:256
        im = im + warp_mat3d_fast(mask, angle2dcm4(0,a,0), 1, depth);
    end
    im(im==0)=NaN;
    
    im = im - minval;
    im = im / maxval;
    imwrite(im*256+1, hot(256), getresultsdir(sprintf('maskanim%.3d.png', counter)));
    counter = counter + 1;
%    imshow(im);colormap(hot);
%     imshow(log(1+im),minmax); colormap(hot);
%     frm = getframe();
%     MOVIE = addframe(MOVIE,frm.cdata);
%   pause(1);
end