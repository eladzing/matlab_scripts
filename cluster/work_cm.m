env

MPSec = 1;
h = 0.7;

RRHOTOT = RHOTOT(MPSec);

%%
[meshY, meshX, meshZ] = meshgrid(1:256, 1:256, 1:256);
meshX = meshX - (256+1)/2;
meshY = meshY - (256+1)/2;
meshZ = meshZ - (256+1)/2;

%%
CMx = sum(meshX(:).*RRHOTOT(:))/sum(RRHOTOT(:))
CMy = sum(meshY(:).*RRHOTOT(:))/sum(RRHOTOT(:))
CMz = sum(meshZ(:).*RRHOTOT(:))/sum(RRHOTOT(:))


%%
cm = [CMx CMy CMz];
new_vr = cart2sphere(Vr(1, cm), 256, cm);
new_rhog = cart2sphere(RHOG(1), 256, cm);