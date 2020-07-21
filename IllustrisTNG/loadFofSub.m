


global DEFAULT_MATFILE_DIR
global simDisplayName
global DRACOFLAG

if ~exist('snap','var')
    fprintf('snap not defined, defaulting to snap=99 \n')
    snap=99;
end

if ~DRACOFLAG
    load([ DEFAULT_MATFILE_DIR '/fofs_subs_' simDisplayName '_snp' num2str(snap) '.mat'])
else
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
end
