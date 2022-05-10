function [subs,fofs,subsInfo]=loadFofSub(snap)
% load the FOF nad subhalo structures either from a local catalog or from
% the simulation directory

global DEFAULT_MATFILE_DIR
global simDisplayName
global DRACOFLAG
global BASEPATH
global MOUNTEDFLAG
if ~exist('snap','var')
    fprintf('snap not defined, defaulting to snap=99 \n')
    snap=99;
end

if ~DRACOFLAG
    try
        load([ DEFAULT_MATFILE_DIR '/fofs_subs_' simDisplayName '_snp' num2str(snap) '.mat'])
    catch
        if MOUNTEDFLAG
            fofs=illustris.groupcat.loadHalos(BASEPATH,snap);
            subs=illustris.groupcat.loadSubhalos(BASEPATH,snap);
        else
            error('%s - could not load catalog. local copy may not exist for snapshot %i. Consider mounting Vera',current_function().upper,snap);
        end
    end
else
    fofs=illustris.groupcat.loadHalos(BASEPATH,snap);
    subs=illustris.groupcat.loadSubhalos(BASEPATH,snap);
end

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);

end