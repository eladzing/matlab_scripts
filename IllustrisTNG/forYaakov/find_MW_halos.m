

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
    
	%HydroCat=illustris.utils.read_catalog('gas','snap',99,'folder','hydrogen');
	readFlag=false;
    
    
end

%%

m200=fofs.Group_M_Crit200.*illUnits.massUnit;
mwMask=log10(m200)>=11.5 & log10(m200)<=12.5;
halInd=find(mwMask);
