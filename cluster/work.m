global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

R = log10(1+RHOG(8));
for a = 0:pi/100:pi/2
    %imshow(warp_mat3d_2(R, angle2dcm4(0,a,0), 1), []); 
    imshow(warp_mat3d_sum(R, angle2dcm4(0,a,0)), []); 
    colormap(hot);
    pause(1/100);
end