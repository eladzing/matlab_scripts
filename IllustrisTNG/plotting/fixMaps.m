%id=[0 1 2 3 4 5 7 9 10 13]; %300
id=[0 1 2 3 4 5 6 8 9 10]; % 100 
%mvs=1.0e+15.*[1.5358  1.3073    1.0333    0.8997    0.8418    0.7335
%0.5812    0.6412    0.6359    0.5539]; % 300
mvs=fofs.Group_M_Crit200(id+1).*illUnits.massUnit; % 100


% open figs
openfig(sprintf('%s/%s',DEFAULT_FIG_DIR,'map_ent_fof10_noVel_TNG100.fig'))
openfig(sprintf('%s/%s',DEFAULT_FIG_DIR,'map_temp_fof10_noVel_TNG100.fig'))
openfig(sprintf('%s/%s',DEFAULT_FIG_DIR,'map_nDens_fof10_noVel_TNG100.fig'))

mv=mvs(id==10)

for i=1:9
    figure(i)
    
    titlemine(mk_mvir_string(mv,'200,c'))
    
    switch i
        case {1, 2, 3}
            cmap=brewermap(256,'*PuOr');
            colormap(cmap)
             caxis([2 3.5])
        case {4,5,6}
            caxis([6.5 8])
        case{ 7, 8, 9}
            cmap=brewermap(256,'YlOrRd');
            colormap(cmap)
            caxis([-5.5 -2.5])
    end
    
end


fof='10';
baseName='map%s_%s_fof%s_TNG100';

printout_fig(1,sprintf(baseName,'YZ','ent',fof))
printout_fig(2,sprintf(baseName,'XZ','ent',fof))
printout_fig(3,sprintf(baseName,'XY','ent',fof))

printout_fig(5,sprintf(baseName,'XZ','temp',fof))
printout_fig(4,sprintf(baseName,'YZ','temp',fof))
printout_fig(6,sprintf(baseName,'XY','temp',fof))

printout_fig(7,sprintf(baseName,'YZ','nDens',fof))
printout_fig(8,sprintf(baseName,'XZ','nDens',fof))
printout_fig(9,sprintf(baseName,'XY','nDens',fof))

close all
