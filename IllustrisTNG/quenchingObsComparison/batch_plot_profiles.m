%% batch plots 


%% radius 

plotRadProf(hprofs50,'CGMall','TNG50',1,'print')
plotRadProf(hprofs100,'Gal','TNG100',1,'print')
plotRadProf(hprofs100,'CGMall','TNG100',1,'print')
close all


%ssfr 
plotsSFRProf(hprofs50,'Gal','TNG50',1,'print')
plotsSFRProf(hprofs50,'CGMall','TNG50',1,'print')
plotsSFRProf(hprofs100,'Gal','TNG100',1,'print')
plotsSFRProf(hprofs100,'CGMall','TNG100',1,'print')
close all

plotRadProf_massN(hprofs50,'Gal','TNG50',1,'print')
plotRadProf_massN(hprofs50,'CGMall','TNG50',1,'print')