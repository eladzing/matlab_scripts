%% reproduce Zhang+2026
%% perliminaries

bp=illustris.set_env(100);
snap=99; 

illustris.utils.set_illUnits(snap);

global illUnits;
global DEFAULT_MATFILE_DIR
global simDisplayName
global cosmoStruct

%%
if ~exist('readFlag','var')
    readFlag=true;
end

% read FOF and SUBFIND data, as well as free-fall time profiles
if readFlag
    fprintf(' *** Reading data *** \n');
    [subs, fofs, subsInfo]=illustris.loadFofSub(snap);

    readFlag=false;
%    fofs = illustris.utils.addTvirFofs( fofs);

end
R200c=double(fofs.Group_R_Crit200(subsInfo.hostFof+1));
%%
% basic mask 
massThresh=5e9; % threshold for *stellar* mass
massAllGals=illustris.utils.get_stellar_mass(subs,'gal');
maskBase=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');

papMask=massAllGals>massThresh & subsInfo.isCentral & subs.SubhaloBHMass>0;% & subs.SubhaloMassType(1,:)>0;

fprintf('Number of gals in base mask : %i \n',sum(maskBase & subs.SubhaloBHMass>0))
fprintf('Number of gals in paper mask : %i \n',sum(papMask))
%% go over galaxies 

indxBase=find(papMask);

ids=indxBase-1;
cnt=0;

step=5;
stepNext=5;

for id=ids
    % for following progress
    perCent=floor((cnt+1)/length(ids)*100);

    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end

    cnt=cnt+1;
    gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas');
    if gas.count==0
        continue
    end

    sfMask=gas.StarFormationRate>0;
    %set mass
    mass=double(gas.Masses.*illUnits.massUnit); %in Solarmass
    gas=illustris.utils.addTemperature(gas);
    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
    gasDist=sqrt(sum(double(gas.newCoord).^2,1));

    %  radial mask 
    rv=double(R200c(id+1));
    
    radMask(1,:)=gasDist<=0.03*rv;
    radMask(2,:)=gasDist<=0.15*rv;
    radMask(3,:)=gasDist<=rv;
    radMask(4,:)=gasDist>rv;

    coldMask=gas.Temperature<1e4 & ~sfMask;
    coolMask=gas.Temperature>=1e4 & gas.Temperature<1e5 & ~sfMask;
    warmMask=gas.Temperature>=1e5 & gas.Temperature<1e6 & ~sfMask;
    hotMask=gas.Temperature>=1e6 & ~sfMask;

        
    for k=1:4
        gasMass.cold(k,cnt)=sum(mass(coldMask & radMask(k,:)));
        gasMass.cool(k,cnt)=sum(mass(coolMask & radMask(k,:)));
        gasMass.warm(k,cnt)=sum(mass(warmMask & radMask(k,:)));
        gasMass.hot(k,cnt)=sum(mass(hotMask & radMask(k,:)));
        gasMass.tot(k,cnt)=sum(mass(radMask(k,:)));
    end
        
clear radMask 
end

%% save to catalog file
global DRACOFLAG
if DRACOFLAG

    fname=sprintf('gasMasses_zhang2026_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'PropStruct','massHistStruct','-v7.3')

    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);

end
fprintf(' *** DONE!  *** \n');