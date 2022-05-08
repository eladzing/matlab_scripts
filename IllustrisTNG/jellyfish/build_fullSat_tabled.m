%% Need build an objectTable similar to the CJF one but which includes all
% satellites which fit the CJF project criteria

snap100=[33 40 50 59 67 72 78 84 91 99];
snap50= [33 40 50 59 67:99];


%% Initialize  table
fullSatTable=table([],[],[],[],...
    'variableNames',{'sim','snap','subfind','tag','hostID'});

%% TNG100
bp=illustris.set_env(100,mountFlag);
massThresh=10^9.5;

indStart=1;
for i=1:length(snap100)
    snap=snap100(i);
    illustris.utils.set_illUnits(snap);
    
    %% load catalogs
    [subs,fofs,subsInfo]=illustris.loadFofSub(snap);
    
    %% generate mask
    smass=illustris.utils.get_stellar_mass(subs);
    
    msk=smass>=massThresh & ~subsInfo.isCentral & subs.SubhaloFlag;
    
    % check that the mask encompass the CJF sample
    sampleMask=objectTable.sim=="TNG100" & objectTable.snap==snap;
    sampleID=objectTable.subfind(sampleMask);
    maskID=find(msk)-1;
    missingInds=[];
    for j=1:length(sampleID)
        if ~any(sampleID(j)==maskID)
            warning('%s - missing index from sample(TNG100, %i) : %i',...
                current_function().upper,snap,sampleID(j));
            missingInds(end+1)=sampleID(j);
        end
    end
    
    %% fill up table
    indEnd=indStart+sum(msk)-1;
    tabInds=indStart:indEnd;
    
    galInds=maskID+1; % indices in the catalogs
    hostInds=subsInfo.hostFof(galInds)+1;
    
    
    fullSatTable.sim(tabInds)="TNG100";
    fullSatTable.snap(tabInds)=snap;
    fullSatTable.subfind(tabInds)=maskID;
    fullSatTable.hostID(tabInds)=subsInfo.hostFof(msk);
    
    %     % build tag
    %     snp=string(['0' num2str(snap)]);
    %     tg1(1:length(mskInds),1)="snp";
    %     tg2(1:length(mskInds),1)="subid";
    %     tg3(1:length(mskInds),1)="typ:";
    %     clsTab.tag=join(['TNG100' tg1 snp tg2 sid tg3 imgType'],'');
    %     fullSatTable.tag(tabInds)=
    
    %% fill up properties
    
    
    
    fullSatTable.hostM200c(tabInds)=illUnits.massUnit.*...
        double(fofs.Group_M_Crit200(hostInds));
    fullSatTable.hostR200c(tabInds)=illUnits.lengthUnit.*...
        double(fofs.Group_R_Crit200(hostInds));
    fullSatTable.hostM200m(tabInds)=illUnits.massUnit.*...
        double(fofs.Group_M_Mean200(hostInds));
    fullSatTable.hostR200m(tabInds)=illUnits.lengthUnit.*...
        double(fofs.Group_R_Mean200(hostInds));
    fullSatTable.stellarMass(tabInds)=illUnits.massUnit.*...
        double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,galInds));
    fullSatTable.galStellarMass(tabInds)=;
    fullSatTable.gasMass(tabInds)=;  % total
    fullSatTable.galGasMass(tabInds)=; % total
    fullSatTable.galSFR(tabInds)=;
    fullSatTable.pos=zeros(3,ngal);
    fullSatTable.vel=zeros(3,ngal);
    fullSatTable.rpos(tabInds)=;
    fullSatTable.vrad(tabInds)=;
    fullSatTable.galBHMass(tabInds)=;
    
    indStart=indEnd+1;
end
