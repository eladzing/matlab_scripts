%% batch plots 


%% radius 
plotRadProf(hprofs50,'Gal','TNG50',1,'xl',[0 2.1],'print')
plotRadProf(hprofs50,'CGMall','TNG50',1,'xl',[0 2.1],'print')
plotRadProf(hprofs50,'Gal','TNG50',3,'xl',[0 2.1],'print')
plotRadProf(hprofs50,'CGMall','TNG50',3,'xl',[0 2.1],'print')

plotRadProf(hprofs100,'Gal','TNG100',1,'xl',[0 2.1])
plotRadProf(hprofs100,'CGMall','TNG100',1,'xl',[0 2.1])
close all


%% ssfr 
plotsSFRProf(hprofs50,'Gal','TNG50',1)
plotsSFRProf(hprofs50,'CGMall','TNG50',1)
plotsSFRProf(hprofs100,'Gal','TNG100',1)
plotsSFRProf(hprofs100,'CGMall','TNG100',1)
close all

plotRadProf_massN(hprofs50,'Gal','TNG50',1)
plotRadProf_massN(hprofs50,'CGMall','TNG50',1)