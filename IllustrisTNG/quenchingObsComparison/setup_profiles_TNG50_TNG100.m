binEdges=0.1:0.2:2.1;

%% load important things TNG50
bp=illustris.set_env(50,'nomount');
snap=99;
[subs50,fofs50,subsInfo50]=illustris.loadFofSub(snap);

massThresh=10^8.3;
satMask50=illustris.infrastructure.generateMask('subs',subs50','fofs',fofs50,'mass',massThresh,'snap',snap,'sats');
%% load HI H2 values by components and set profiles (3d)

global DEFAULT_MATFILE_DIR
global simDisplayName
load([DEFAULT_MATFILE_DIR '\hih2Catalog_snp99_' simDisplayName '.mat'])
hih2Struct50=hih2Struct;

hprofs50=generate_hih2_profilesFull_3d(hih2Struct50,fofs50,subs50,'mask',satMask50,'mean200','radBins',binEdges);


clear hih2Struct
%% load important things TNG100

bp=illustris.set_env(100,'nomount');
snap=99;
[subs100,fofs100,subsInfo100]=illustris.loadFofSub(snap);

massThresh=10^9;
satMask100=illustris.infrastructure.generateMask('subs',subs100','fofs',fofs100,'mass',massThresh,'snap',snap,'sats');
%% load HI H2 values by components and set profiles (3d)


load([DEFAULT_MATFILE_DIR '\hih2Catalog_snp99_' simDisplayName '.mat'])
hih2Struct100=hih2Struct;

hprofs100=generate_hih2_profilesFull_3d(hih2Struct100,fofs100,subs100,'mask',satMask100,'mean200','radbins',binEdges);
