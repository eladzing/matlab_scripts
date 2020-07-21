%%
work_stream2
SS = S(8);
TT= T(8);
new_verts = filter_verts(verts, 40, SS, 0.1,0.5);
new_verts2 = filter_verts(new_verts, 40, log10(TT), 0.5,inf);
mask = f(new_verts2,40);
idxs = find(mask);

%%
clear mask SS TT meshZ Vrr Vxx Vyy Vzz cosa meshX meshY Vabs

%%
MPSec = 8;
[result_map result_vrho result_vtemp result_count] = work_bird(MPSec, idxs);

%%
figure;
subplot(2,2,2);
work_bird_plots;
[xv, yv] = getline(gcf);
line(xv,yv, 'Color','g');

rho_log = log10(RHOG(MPSec));
temp_log = log10(T(MPSec));
mask = inpolygon(temp_log,rho_log,xv,yv);

plot_bird2slice(MPSec,0,0,mask);

print('-dpng','-r400', getresultsdir(sprintf('stream %s.png',datestr(now,'yyyy-mm-dd HH-MM-SS'))))
