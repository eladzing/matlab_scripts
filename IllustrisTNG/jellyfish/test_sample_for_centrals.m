

global illUnits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat'])
snaps=unique(objectTable.snap);
sims=unique(objectTable.sim);
%% initialize

% snap & sim and id
centralTest.snap=objectTable.snap;
centralTest.sim=objectTable.sim;
centralTest.subfind=objectTable.subfind;
centralTest.tag=objectTable.tag;

ngal=length(objectTable.subfind);

centralTest.isCentral=false(1,ngal);
centralTest.rad10=false(1,ngal);
centralTest.rad5=false(1,ngal);
centralTest.done=false(1,ngal);
%%run over simulations and snapshots
for k=1:length(sims)
    
    bp=illustris.set_env(sims(k));

    global LBox
    fprintf('running on simulation %s \n ', sims(k));
    
    fprintf('going over %i snapshots \n',length(snaps));
    simMask=strcmp(objectTable.sim,sims(k));
    
    for i=1:length(snaps)
        
        snap=snaps(i);
        fprintf('%i )  snap %i ',i,snap);
        % find indices of gals from the snapshot
        snapMask=objectTable.snap==snap;
        %fullMask=snapMask & simMask;
        inds=find(snapMask & simMask);
        
        if isempty(inds) % skip snapshot with no relevant objects
            fprintf(' .... no objects \n');
            continue
        else
            fprintf(' ... %i objects \n',length(inds));
        end
        
        fprintf('   Loading catalogs of snap %i \n',snap);
        [subs,fofs,subsInfo]=illustris.loadFofSub(snap);
        illustris.utils.set_illUnits(snap);
        
        fprintf('   Getting galaxy properties \n');
        
        
        
        galInds=objectTable.subfind(inds)+1; % indices in the catalogs
        hostInds=subsInfo.hostFof(galInds)+1;% indices in the catalogs
        
        %% fill out the properties
        centralTest.done(inds)=true;
        centralTest.isCentral(inds)=subsInfo.isCentral(galInds);
        
        
        
         %% position data
        galPos=double(subs.SubhaloPos(:,galInds)); % global position in simulation box, in simulation units
        hostPos=double(fofs.GroupPos(:,hostInds));% global position of host in simulation box, in simulation units
        
        rpos=findDistance(galPos,hostPos,LBox,3);
        r200=illUnits.lengthUnit.*double(fofs.Group_R_Crit200(hostInds));
        
        rp=rpos./r200;
        
        centralTest.rad10(inds)=rp<=0.1;
        centralTest.rad5(inds)=rp<=0.5;
        
    end
end


fname=sprintf('jf_testCentral_CJF.mat');
save([DEFAULT_MATFILE_DIR '/' fname],'centralTest','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);





% centralTest.vel(inds,3)=...
%
%     centralTest.rpos(inds)=...
%
%     centralTest.vrad(inds)=...
%
%     centralTest.vtan(inds)=...

%massThresh=10^8.5; % threshold for *stellar* mass




