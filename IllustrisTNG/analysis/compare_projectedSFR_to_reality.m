%% start at redshift 0 
snap=99;
% load galaxy data 
   global DRACOFLAG
        if DRACOFLAG
            fofs=illustris.groupcat.loadHalos(bp,snap);
            subs=illustris.groupcat.loadSubhalos(bp,snap);
    
        else
            loadFofSubTNG100
        end

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
illustris.utils.set_illUnits(snap)
global illUnits

% identify relavent central galaxies and set mask 
massThresh=10^9; % threshold for *stellar* mass
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
spbMask = mk_splashback_mask('time',5,'both',0.1);

mask=galMask & ~spbMask; 
%centralMask= subsInfo.isCentral(tCoolStruct.galMask);

% load trees and follow main branch to snap 91. 

global BASEPATH
treeType='SubLink_gal';
treeFields={'SnapNum','SubfindID','SubhaloVel','SubhaloPos','SubhaloHalfmassRadType',...
    'SubhaloID','SubhaloMassInRadType','SubhaloSFRinRad','Group_M_Crit200','SubfindID',...
    'SubhaloHalfmassRadType','FirstSubhaloInFOFGroupID','GroupPos','Group_R_Crit200'};

ids=find(mask);

for i=inds

tree = illustris.sublink.loadTree(BASEPATH,snap,i-1,treeFields,true,treeType);


% match to snap 99 galaxies 