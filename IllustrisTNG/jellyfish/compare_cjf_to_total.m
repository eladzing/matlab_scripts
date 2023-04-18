%% TNG50 
bp=illustris.set_env(50);
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '\cosmic_jellyfish_objectTable.mat']);
load([DEFAULT_MATFILE_DIR '/jf_galProperties_CJF.mat']);

snaps=unique(objectTable.snap(objectTable.sim=="TNG50"));

for k=1:length(snaps) 
   [subs,fofs,subsInfo]=illustris.loadFofSub(snaps(k));
   ms=illustris.utils.get_stellar_mass(subs); % find stellar mass 
   mask=ms>=10^8.3 & ~subsInfo.isCentral; % find sat no. 
   hostid=unique(subsInfo.hostFof(mask)); % hosts of satellites
   
   outStruct.snap50(k)=snaps(k);
   outStruct.totalSatNum50(k)=sum(mask);
   outStruct.totalHostNum50(k)=length(hostid);
   
   
   mskCJF=objectTable.sim=="TNG50" & objectTable.snap==snaps(k);
     
   outStruct.cjfSatNum50(k)=sum(mask);
   outStruct.cjfhostNum50(k)=length(unique(galProps.hostID(msk)));
end

%% TNG100 
bp=illustris.set_env(100);

snaps=unique(objectTable.snap(objectTable.sim=="TNG100"));

for k=1:length(snaps) 
   [subs,fofs,subsInfo]=illustris.loadFofSub(snaps(k));
   ms=illustris.utils.get_stellar_mass(subs); % find stellar mass 
   mask=ms>=10^9.5 & ~subsInfo.isCentral; % find sat no. 
   hostid=unique(subsInfo.hostFof(mask)); % hosts of satellites
   
   outStruct.snap100(k)=snaps(k);
   outStruct.totalSatNum100(k)=sum(mask);
   outStruct.totalHostNum100(k)=length(hostid);
   
   
   mskCJF=objectTable.sim=="TNG100" & objectTable.snap==snaps(k);
     
   outStruct.cjfSatNum100(k)=sum(mask);
   outStruct.cjfhostNum100(k)=length(unique(galProps.hostID(msk)));
end