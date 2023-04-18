%% TNG50 
bp=illustris.set_env(50);
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat']);
load([DEFAULT_MATFILE_DIR '/jf_galProperties_CJF.mat']);

snaps=unique(objectTable.snap(objectTable.sim=="TNG50"));
fprintf('snap = ');
for k=1:length(snaps)
    fprintf('%i, ',snaps(k));
   [subs,fofs,subsInfo]=illustris.loadFofSub(snaps(k));
   ms=illustris.utils.get_stellar_mass(subs); % find stellar mass 
   mask=ms>=10^8.3 & ~subsInfo.isCentral; % find sat no. 
   hostid=unique(subsInfo.hostFof(mask)); % hosts of satellites
   
   outStruct.snap50(k)=snaps(k);
   outStruct.totalSatNum50(k)=sum(mask);
   outStruct.totalHostNum50(k)=length(hostid);
   
   
   mskCJF=objectTable.sim=="TNG50" & objectTable.snap==snaps(k);
     
   outStruct.cjfSatNum50(k)=sum(mskCJF);
   outStruct.cjfhostNum50(k)=length(unique(galProps.hostID(mskCJF)));
end
fprintf(';\n');
%% TNG100 
bp=illustris.set_env(100);

snaps=unique(objectTable.snap(objectTable.sim=="TNG100"));
fprintf('snap = ');
for k=1:length(snaps) 
    fprintf('%i, ',snaps(k));
   [subs,fofs,subsInfo]=illustris.loadFofSub(snaps(k));
   ms=illustris.utils.get_stellar_mass(subs); % find stellar mass 
   mask=ms>=10^9.5 & ~subsInfo.isCentral; % find sat no. 
   hostid=unique(subsInfo.hostFof(mask)); % hosts of satellites
   
   outStruct.snap100(k)=snaps(k);
   outStruct.totalSatNum100(k)=sum(mask);
   outStruct.totalHostNum100(k)=length(hostid);
   
   
   mskCJF=objectTable.sim=="TNG100" & objectTable.snap==snaps(k);
     
   outStruct.cjfSatNum100(k)=sum(mskCJF);
   outStruct.cjfhostNum100(k)=length(unique(galProps.hostID(mskCJF)));
end
fprintf(';\n');
save([DEFAULT_MATFILE_DIR '/cjf_total_sample_comp.mat'],'outStruct');
