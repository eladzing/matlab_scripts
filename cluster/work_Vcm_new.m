%% Set environment
env

%% 18 Boxes, or 18*2 half boxes
Vxx = cart2sphere(Vx(1), 18*2);
Vyy = cart2sphere(Vy(1), 18*2);
Vzz = cart2sphere(Vz(1), 18*2);
RRhog = cart2sphere(RHOG(1), 18*2);
ds = ds_sphere(1);

Vxx = Vxx(1:36,:,:);
Vyy = Vyy(1:36,:,:);
Vzz = Vzz(1:36,:,:);
ds = ds(1:36,:,:);
RRhog = RRhog(1:36,:,:);

RRhog = RRhog.*ds;
M = sum(RRhog(:));

VcmX = sum(Vxx(:).*RRhog(:))/M
VcmY = sum(Vyy(:).*RRhog(:))/M
VcmZ = sum(Vzz(:).*RRhog(:))/M

%save('Vcm', 'VcmX', 'VcmY', 'VcmZ');