%assumes mask (256^3) is defined

%%
global xmesh;
global ymesh;
global zmesh;
[xmesh ymesh zmesh] = meshgrid(1:256,1:256,1:256);

%%

% MOVIE = avifile('anim.avi','fps',2);%,'COMPRESSION','Cinepak');

%[min(log(1+im(:))) max(log(1+im(:)))]

% %for masked gas density - use LOG10 (!)
% load verygoodmask
% original = RHOG(8);
% original(~mask) = 0;
% minval = 10.35;
% maxval = 14;

% %for dark matter - use LOG10 (!)
% original = RHODM(8);
% minval = 9;%7.15;
% maxval = 15.5;

% %for VR - - dont use LOG10 (!)
% load verygoodmask
% original = Vr(8);
% original(~mask)=0;
% %original(original<0)=0; %to see the outliars in speed
% original(original>0)=0;
% original = -original;
% minval = 0;
% maxval = 50000;

% %for RHOG in a=0.6
% original = RHOG(8);
% minval = 12.2890;
% maxval = 15;

%for VR - - dont use LOG10 (!)
load verygoodmask
original = Vr(8).*RHOG(8);
original(~mask)=0;
%original(original<0)=0; %to see the outliars in speed
original(original>0)=0;
original = -original;
minval =9;
maxval = 18;

counter = 0;
for a = 0:pi/100:2*pi%0:pi/5:pi%0:pi/2:2*pi
    a/pi*180
    im = zeros(256,256);
    for depth = 1:1:256
        im = im + warp_mat3d_fast(original, angle2dcm4(0,a,0), 1, depth);
    end

    im = log10(im);
    im(isinf(im)) = NaN;
    
%        minval = min(im(:))
%        maxval = max(im(:))
    
    im = im - minval;
    im = im / (maxval-minval);
    imwrite(im*256, hot(256), getresultsdir(sprintf('flux%.3d.png', counter)));
    counter = counter + 1;
%    imshow(im);colormap(hot);
%     imshow(log(1+im),minmax); colormap(hot);
%     frm = getframe();
%     MOVIE = addframe(MOVIE,frm.cdata);
%   pause(1);
end