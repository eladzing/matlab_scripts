

global illUnits
global DEFAULT_MATFILE_DIR

load([DEFAULT_MATFILE_DIR '/jf_objectTable_TNG50.mat']);
snaps=unique(objectTable.snap);

%% initialize
ngal=length(objectTable.subfind);
galProps.stellarMass=zeros(1,ngal);
galProps.hostM200c=zeros(1,ngal);
galProps.hostR200c=zeros(1,ngal);
galProps.hostM200m=zeros(1,ngal);
galProps.hostR200m=zeros(1,ngal);
galProps.galStellarMass=zeros(1,ngal);
galProps.gasMass=zeros(1,ngal);
galProps.galGasMass=zeros(1,ngal);
galProps.galSFR=zeros(1,ngal);
galProps.pos=zeros(3,ngal);
galProps.vel=zeros(3,ngal);
galProps.rpos=zeros(1,ngal);
galProps.vrad=zeros(1,ngal);
galProps.vtan=zeros(1,ngal);
%galProps.nSatHost=zeros(1,ngal);
%galProps.rhoLocal=zeros(size(objectTable.subfind));
%galProps.ramPress=zeros(size(objectTable.subfind));
%galProps.mach=zeros(size(objectTable.subfind));

for k=1:length(sims)
    
    
    for i=1:length(snaps)
    
    snap=snaps(i);
    
    fprintf('Loading catalogs of snap %i \n',snap);
    loadFofSub;
    illustris.utils.set_illUnits(snap);
    
    
    fprintf('Getting galaxy properties \n');
    
    % find indices of gals from the snapshot 
    inds=find(objectTable.snap==snap);
    
    galIDs=objectTable.subfind(inds)+1;
    hostIDs=objectTable.host(inds)+1;
    
    %% fill out the properties 
    galProps.stellarMass(inds)=illUnits.massUnit.*...
        subs.SubhaloMassType(illustris.partTypeNum('stars')+1,galIDs);

    galProps.hostM200c(inds)=illUnits.massUnit.*...
        fofs.Group_M_Crit200(hostIDs);
    
    galProps.hostR200c(inds)=illUnits.lengthUnit.*...
        fofs.Group_R_Crit200(hostIDs);
    
    galProps.galStellarMass(inds)=illUnits.massUnit.*...
        subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galIDs);
   
    galProps.gasMass(inds)=illUnits.massUnit.*...
        subs.SubhaloMassType(illustris.partTypeNum('gas')+1,galIDs);

    galProps.galGasMass(inds)=illUnits.massUnit.*...
        subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,galIDs);

    galProps.galBHMass(inds)=illUnits.massUnit.*...
        subs.SubhaloMassInRadType(illustris.partTypeNum('BH')+1,galIDs);
    
     galProps.galSFR(inds)=...
        subs.SubhaloSFRinRad(galIDs);
    
    % More to come... 
    
end

fname=sprintf('jf_galProperties_%s.mat',simDisplayName);
save([DEFAULT_MATFILE_DIR '/' fname],'galProps','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);





% galProps.vel(inds,3)=...
%     
%     galProps.rpos(inds)=...
% 
%     galProps.vrad(inds)=...
%     
%     galProps.vtan(inds)=...

%massThresh=10^8.5; % threshold for *stellar* mass




