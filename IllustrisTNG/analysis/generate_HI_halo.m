%% This Script generates HI and H_2 masses for eash subhalo - also breaking it down to CGM components  

%% set framework
% the simulation, snapshot, and environment should be set prior to running
% the script

illustris.utils.set_illUnits(snap);

global illUnits
global DEFAULT_MATFILE_DIR
global simDisplayName

if ~exist('readFlag','var')
    readFlag=true;
end

% read FOF and SUBFIND data, as well as free-fall time profiles
if readFlag
    fprintf(' *** Reading data *** \n');
    
	fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
	HydroCat=illustris.utils.read_catalog('gas','snap',99,'folder','hydrogen');
	readFlag=false;
    
    
end

% load a given halo 
id = ;
fields={};
gas=illustris.snapshot.loadSubhalo(bp,snap,id,'gas',fields);

% get the relevant particle indices
subset = illustris.snapshot.getSnapOffsets(bp,snap,id,'Subhalo');

%particle offests are: 
firstInd=subset.offsetType(1)+1;
lastInd=firstInd+int64(subset.lenType(1))-1;

% get the total neutral mass
gas.mH_type={'BR','GK','KMT'};
gas.mh=HydroCat.MH(firstInd:lastInd)';

gas.mH2(1,:)=HydroCat.MH2BR(firstInd:lastInd)';
gas.mH2(2,:)=HydroCat.MH2GK(firstInd:lastInd)';
gas.mH2(3,:)=HydroCat.MH2KMT(firstInd:lastInd)';
gas.mHi(1,:)=gas.mh-HydroCat.MH2BR(firstInd:lastInd)';
gas.mHi(2,:)=gas.mh-HydroCat.MH2GK(firstInd:lastInd)';
gas.mHi(3,:)=gas.mh-HydroCat.MH2KMT(firstInd:lastInd)';




