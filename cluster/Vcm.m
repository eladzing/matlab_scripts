function result = Vcm(MPSec)

Vxx = Vx(MPSec);
Vyy = Vy(MPSec);
Vzz = Vz(MPSec);
Rhog = RHOG(MPSec);
M = sum(Rhog(:));

result = [sum(Vxx(:).*Rhog(:)) sum(Vyy(:).*Rhog(:)) sum(Vzz(:).*Rhog(:))]/M;