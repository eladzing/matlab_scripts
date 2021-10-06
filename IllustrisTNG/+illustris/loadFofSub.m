function [subs,fofs,subsInfo]=loadFofSub(snap)
% load the FOF nad subhalo structures either from a local catalog or from
% the simulation directory

global DEFAULT_MATFILE_DIR
global simDisplayName
global DRACOFLAG
global BASEPATH

if ~exist('snap','var')
    fprintf('snap not defined, defaulting to snap=99 \n')
    snap=99;
end

if ~DRACOFLAG
    load([ DEFAULT_MATFILE_DIR '/fofs_subs_' simDisplayName '_snp' num2str(snap) '.mat'])
else
    fofs=illustris.groupcat.loadHalos(BASEPATH,snap);
    subs=illustris.groupcat.loadSubhalos(BASEPATH,snap);
end

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);

end