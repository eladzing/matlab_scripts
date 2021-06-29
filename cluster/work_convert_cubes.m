global FILE_FORMAT;
FILE_FORMAT = '../data/%s_a1.001L%dMpc.dat';

% res = cart2sphere(Vr(4));
% save(sprintf('../data/Vr_sphere_%d', 4), 'res');
% 
% for MPSec = [8]
%     res = cart2sphere(RHOG(MPSec));
%     save(sprintf('../data/RHOG_sphere_%d', MPSec), 'res');
%     res = cart2sphere(RHOTOT(MPSec));
%     save(sprintf('../data/RHOTOT_sphere_%d', MPSec), 'res');
%     res = cart2sphere(T(MPSec));
%     save(sprintf('../data/T_sphere_%d', MPSec), 'res');
%     res = cart2sphere(S(MPSec));
%     save(sprintf('../data/S_sphere_%d', MPSec), 'res');
%     res = cart2sphere(Vr(MPSec));
%     save(sprintf('../data/Vr_sphere_%d', MPSec), 'res');
% end

%%
for MPSec = [8 4 2 1]
    res = cart2sphere(Vr(MPSec));
    save(sprintf('../data/Vr_sphere_%d', MPSec), 'res');
end