% env
% MPSec = 1;
% [vx, vy, vz, v1,v2,v3] = V_Vcm(MPSec);
% 
% small_len = length(vx);
% grid_vec = (1:small_len)-((small_len+1)/2);
% [xxx,yyy,zzz] = meshgrid(grid_vec,grid_vec,grid_vec);
% radii = sqrt(xxx.^2+yyy.^2+zzz.^2);
% clear xxx yyy zzz
% 
% units

% for radius = [25 50]
%     radius
%     mask = (radii > radius);
% %     convn(vx, kernel, 'same')
% %     convn(vy, kernel, 'same')
% %     convn(vz, kernel, 'same')
%     vx2 = ifftn(ifftshift(fftshift(fftn(vx)).*double(mask)));
%     vy2 = ifftn(ifftshift(fftshift(fftn(vy)).*double(mask)));
%     vz2 = ifftn(ifftshift(fftshift(fftn(vz)).*double(mask)));
% 
%     new_v_sq = abs(vx2).^2 + abs(vy2).^2 + abs(vz2).^2;
% 
% 
%     %this line calculates the pressure locally. instead we can calculate only
%     %the V_sq profile and multiply it by RHOG profile
%     %pressure = 0.5*new_v_sq.*1e3.*RHOG(1)*(Msun/((MPc)^3)); %in MKS!
%     pressure = 0.5*new_v_sq.*1e6;
%     pressure_sphere = cart2sphere(pressure);
%     pressure_profile{radius}  = squeeze(sum(sum(pressure_sphere,3),2))'/(256^2);
% end

% mask = (radii < 10);
% [xx yy zz] = ind2sub(size(mask),find(mask));
% %[yyy,xxx,zzz] = meshgrid(grid_vec,grid_vec,grid_vec);
% [yyy,xxx,zzz] = meshgrid(-128:127,-128:127,-128:127);
% th=10;
% xxx = xxx(128-th:128+th,128-th:128+th,128-th:128+th);
% yyy = yyy(128-th:128+th,128-th:128+th,128-th:128+th);
% zzz = zzz(128-th:128+th,128-th:128+th,128-th:128+th);
% 
% vx2 = zeros(size(xxx));
% vy2 = zeros(size(xxx));
% vz2 = zeros(size(xxx));
% len_ball = length(xx(:));
% tic
% for idx = 1:length(xxx(:))
%     xx2 = xx - xxx(idx);
%     yy2 = yy - yyy(idx);
%     zz2 = zz - zzz(idx);
%     %idxs = sub2ind(size(mask), xx,yy,zz);
%     idxs = xx2 + (yy2-1)*256 + (zz2-1)*(256^2);
% 
%     vx2(idx) = mean(vx(idxs));
%     vy2(idx) = mean(vy(idxs));
%     vz2(idx) = mean(vz(idxs));
% end
% toc