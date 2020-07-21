function set_illUnits(snap) 
% SET_ILLUNITS set the conversion factors from TNG simulation units to more
%reasonable things, also changes from comoving to physical units. 

global illUnits 


% if no snap given, use z=0;
if ~exist('snap','var')
    snap=99;
end

illUnits=illustris.utils.calc_illUnits(snap);  
    
    
end