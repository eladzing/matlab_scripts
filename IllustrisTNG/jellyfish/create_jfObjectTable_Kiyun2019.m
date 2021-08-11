%% short script to generate a jellyfish object table from Kiyun's catalogs 

snaps=[63 72 84 99];
len=0;

subfindID=[];
snapArr=[];
scoreArr=[];
for i=1:length(snaps)
    
    snap=snaps(i);
    
    scores=illustris.utils.read_Kiyun_jellyfish_catalog(snap);
    msk=scores>=0;
    
    subID=(find(msk)-1)';
    snp=ones(sum(msk),1).*snap;
    scr=scores(msk)';
    
    subfindID=cat(1,subfindID,subID);
    snapArr=cat(1,snapArr,snp);
    scoreArr=cat(1,scoreArr,scr);
    
    jfStruct(i).subID=subID;
    jfStruct(i).snap=snp;
    jfStruct(i).scores=scr;
    len=len+sum(msk);
end
 
%% build objectTable 
objectTable=table(subfindID,snapArr,scoreArr,...
    'variableNames',{'subfind','snap','score'});
global DEFAULT_MATFILE_DIR
save([DEFAULT_MATFILE_DIR '/yun19_jellyfish_objectTable.mat'],'objectTable','jfStruct');