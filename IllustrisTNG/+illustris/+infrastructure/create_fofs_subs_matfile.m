snaps=[33    40    50    67    84    99];

global DEFAULT_MATFILE_DIR
global simDisplayName

for i=1:length(snaps)
    
    snap=snaps(i);
    fprintf('snap %i \n',snap);
    fprintf('reading fofs \n');
    
    fofs=illustris.groupcat.loadHalos(bp,snap);
    fprintf('reading subs \n');
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
    fname=['fofs_subs_' simDisplayName '_snp' num2str(snap) '.mat'];
    fprintf('writing to file: %s  \n',fname);
    save([DEFAULT_MATFILE_DIR '/' fname],'fofs','subs', '-v7.3')
end