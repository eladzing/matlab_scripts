%% load HI-H2 catalog 
snap=99;
HydroCat=illustris.utils.read_catalog('gas','snap',99,'folder','hydrogen');

% load a given halo 
id = 554728;
fields={};
gas=illustris.snapshot.loadSubhalo(bp,snap,id,'gas',fields);

% get the relevant particle indices
subset = illustris.snapshot.getSnapOffsets(bp,snap,id,'Subhalo');

%particle offests are: 
firstInd=subset.offsetType(1)+1;
lastInd=subset.lenType(1)+1;

% get the total neutral mass
gas.mh=gs.MH(firstInd:lastInd);
gas.mHi_BR=mh-gs.MH2BR(firstInd:lastInd);
gas.mHi_GK=mh-gs.MH2GK(firstInd:lastInd);
gas.mHi_KMT=mh-gs.MH2KMT(firstInd:lastInd);


