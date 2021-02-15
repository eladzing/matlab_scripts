%% prepare list of stellar masses  with central / satellite labels for Joe.
snap=99;

for k=1:3
    switch k
        case 1
            bp=illustris.set_env(100);
        case 2
            bp=illustris.set_env(300);
        case 3
            bp=illustris.set_env(50);
    end
    
    
    global cosmoStruct
    global illUnits
    global simDisplayName
    
    massThresh0=10*[5.7e4 9.4e5 7.6e6]./cosmoStruct.hub;
    ord=10.^(floor(log10(massThresh0)));
    massThresh=floor(massThresh0./ord).*ord;
    
    %% load data
    loadFofSub;
    
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    
    switch simDisplayName
        case 'TNG50'
            massth=massThresh(1);
        case 'TNG100'
            massth=massThresh(2);
        case 'TNG300'
            massth=massThresh(3);
    end
    
    StellarMass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:)...
        .*illUnits.massUnit; % stellar mass within 2*rhalf
   fprintf('generate mask with mass threshold (10^6): %s \n',num2str(massth/1e6));
   mask=illustris.infrastructure.generateMask('subs',subs,'fofs',fofs,...
        'mass',massth);
    
    
    
    outStruct.StellarMass=double(StellarMass(mask));
    outStruct.isCentral=int8(subsInfo.isCentral(mask));
    outStruct.HostMass=double(fofs.Group_M_Crit200(subsInfo.hostFof(mask)+1)).*illUnits.massUnit;
    outStruct.SFR=double(subs.SubhaloSFRinRad(mask));
    
    
    fprintf('total number of galaxies: %s \n',num2str(sum(mask)));
    fprintf('Centrals comprise: %s \n',num2str(sum(outStruct.isCentral)));
    
    %% write to catalog
    catName=['StellarMassCatalog_' simDisplayName];
    folder='StellarMassCatalogs';
    
    
    illustris.utils.write_catalog(outStruct,snap,'name',catName,...
        'path','default','folder',folder,'v');
    
end
