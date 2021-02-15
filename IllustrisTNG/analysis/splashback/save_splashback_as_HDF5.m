%% save splashback catalog as HDF5 

global simDisplayName
global DEFAULT_MATFILE_DIR

snap=99;


%% load splashback catalogs 
fname=sprintf('splashbackCatalog_byRedshift_%s',simDisplayName);
load([DEFAULT_MATFILE_DIR '/' fname]);

fprintf(' *** loading : %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);


fname=sprintf('splashbackCatalog_byTime_%s',simDisplayName);
load([DEFAULT_MATFILE_DIR '/' fname]);
fprintf(' *** loading : %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);


%% generate mask 
outputStruct.done=int8(all(splashbackRedshift.isSat~=-1,1));

%% generate output arrays
outputStruct.RedshiftBins=splashbackRedshift.redshift;
outputStruct.LookbackTimeBins=splashbackTime.Lookback;

% redshift
outputStruct.SubhaloIsSatelliteRedshift=splashbackRedshift.isSat;
outputStruct.SubhaloBeyond05percentR200cRedshift=splashbackRedshift.isFar_5prc;
outputStruct.SubhaloBeyond10percentR200cRedshift=splashbackRedshift.isFar_10prc;
outputStruct.SubhaloBeyond25percentR200cRedshift=splashbackRedshift.isFar_25prc;
outputStruct.SubhaloBeyond50percentR200cRedshift=splashbackRedshift.isFar_50prc;

%time 
outputStruct.SubhaloIsSatelliteTime=splashbackTime.isSat;
outputStruct.SubhaloBeyond05percentR200cTime=splashbackTime.isFar_5prc;
outputStruct.SubhaloBeyond10percentR200cTime=splashbackTime.isFar_10prc;
outputStruct.SubhaloBeyond25percentR200cTime=splashbackTime.isFar_25prc;
outputStruct.SubhaloBeyond50percentR200cTime=splashbackTime.isFar_50prc;

%% save to HDF5 
                
 catName='Subhalos_BackSplash';
       
 folder=['BackSplashFlags/' simDisplayName];
 
 illustris.utils.write_catalog(outputStruct,snap,'name',catName,...
     'path','default','folder',folder,'v');











