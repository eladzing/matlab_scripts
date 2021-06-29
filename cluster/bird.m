function bird(MPSec)

tm=T(MPSec);
ro=RHOG(MPSec);
trop=hist2d(tm,ro,ro,1);

imagesc(2:9,-7:-18,log10(trop(:,:,1)));
title(sprintf('T-ro mass histogram (%dMpc cube)',MPSec))