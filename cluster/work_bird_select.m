%load ('/home/alf/giorae/results/2008-01-10/coldflows/parametric 2008-01-10 16-20-49.mat') %comment this line for interactive selection

load (sprintf('bird_%d',MPSec));

%env06
%global RHOG_RESULT_8
%global T_RESULT_8
%clear RHOG_RESULT_8 T_RESULT_8

figure;
subplot(2,2,2);
work_bird_plots;
[xv, yv] = getline(gcf); % uncomment this line for interactive selection
line(xv,yv, 'Color','g');

rho_log = log10(RHOG(MPSec));
temp_log = log10(T(MPSec));
mask = inpolygon(temp_log,rho_log,xv,yv);

plot_bird2slice(MPSec,thick,0,0,mask);

%thedate = datestr(now,'yyyy-mm-dd HH-MM-SS');
%print('-dpng','-r400', getresultsdir(sprintf('parametric %s.png',thedate)))
%save(getresultsdir(sprintf('parametric %s',thedate)), 'xv', 'yv', 'MPSec', 'thick');
