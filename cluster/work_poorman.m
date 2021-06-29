

for i = 260:-1:10
imshow(squeeze(F8(i,:,:)), [-2 2]*1e11); colormap(jet); colorbar;
title(sprintf('R = %.2f (R/Rvir)', R_Profile(i)/RVIR));
pause(1/10)
end

