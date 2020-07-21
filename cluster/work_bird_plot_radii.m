%velocities are only the sum of all velocities multiplied by density. We
%need to divide by total density to get the actual velocities
velo_x = result_vtemp./result_map;
velo_y = result_vrho./result_map;

velo_norm = sqrt(velo_x.^2+velo_y.^2);
[xx yy] = meshgrid(2:7/99:9,7:11/99:18);
%quiver(xx,yy,velo_x./velo_norm,velo_y./velo_norm,1)
figure;
%velocities are shown in log-scale but with +1, to avoid negative values
quiver(xx,yy,velo_x./velo_norm.*log10(velo_norm+1),velo_y./velo_norm.*log10(velo_norm+1),3, 'm')
hold on
mapcolormap
%result_radii(result_radii==0)=NaN;
imagesc(2:9,7:18,log10(result_radii));colormap(maxcolormap);
%imagesc(2:9,7:18,result_radii, [0 160]);colormap(maxcolormap);
hold on
log_velo_norm = log10(velo_norm);
log_velo_norm(isinf(log_velo_norm)) = NaN;
log_velo_norm = log_velo_norm - min(log_velo_norm(:));
%quiver(xx,yy,velo_x./velo_norm.*log10(velo_norm+1),velo_y./velo_norm.*log10(velo_norm+1),3, 'c')
%good = find(~isnan(velo_x));
%quiver(xx(good),yy(good),velo_x(good)./velo_norm(good).*log_velo_norm(good),velo_y(good)./velo_norm(good).*log_velo_norm(good),1, 'c')
hold on
streamslice(xx,yy,velo_x,velo_y,10)
xlim([2 9]);ylim([7 18]);
title(sprintf('Parametric diagram (%dMpc cube), with average radius as color',MPSec))
xlabel('log T');
ylabel('log \rho_g');
colorbar;

saveas(gcf,getresultsdir(sprintf('parametric with radius (log)%dMpc.png',MPSec)))
