
function work_bird_select_any(res)
%load ('/home/alf/giorae/results/2008-01-10/coldflows/parametric 2008-01-10 16-20-49.mat') %comment this line for interactive selection
global HALO_PATH
load (sprintf('%s/bird%d',HALO_PATH,res));
thick=128;
%env06

%global T_RESULT_8
%clear RHOG_RESULT_8 T_RESULT_8



figure;
subplot(2,2,2);
imagesc(2:9,7:18,log10(result_map));colormap(jet);
%bird(res);
[xv, yv] = getline(gcf); % uncomment this line for interactive selection
line(xv,yv, 'Color','g');

rho_log = log10(RHOG(res));
temp_log = log10(T(res));
mask = inpolygon(temp_log,rho_log,xv,yv);


plot_bird2slice(res,thick,0,0,mask);

%thedate = datestr(now,'yyyy-mm-dd HH-MM-SS');
%print('-dpng','-r400', getresultsdir(sprintf('parametric %s.png',thedate)))
%save(getresultsdir(sprintf('parametric %s',thedate)), 'xv', 'yv', 'res', 'thick');
