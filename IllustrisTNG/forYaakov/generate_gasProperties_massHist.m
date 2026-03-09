%% This Script generates interesting gas properties for the subfind catalogs
% and saves them as catalogs for TNG


%% set framework
% the simulation, snapshot, and environment should be set prior to running
% the script

illustris.utils.set_illUnits(snap);



global illUnits
global DEFAULT_MATFILE_DIR
global simDisplayName
global cosmoStruct


if ~exist('readFlag','var')
    readFlag=true;
end

% read FOF and SUBFIND data, as well as free-fall time profiles
if readFlag
    fprintf(' *** Reading data *** \n');
    [subs, fofs, subsInfo]=illustris.loadFofSub(snap);

    % fofs=illustris.groupcat.loadHalos(bp,snap);
    % subs=illustris.groupcat.loadSubhalos(bp,snap);
    readFlag=false;
    fofs = illustris.utils.addTvirFofs( fofs);
    %load([DEFAULT_MATFILE_DIR '/freeFallTime_profiles_snp' num2str(illUnits.snap) '_' simDisplayName '.mat']);
end

% assign a lower stellar mass threshold based on the TNG box.
if contains(simDisplayName,'100') || contains(simDisplayName,'300')
    massThresh=10^9; % threshold for *stellar* mass
elseif contains(simDisplayName,'50')
    massThresh=10^(8.3); % threshold for *stellar* mass
    massThreshTop=10^(9.5);
else
    error('mass threshold not set  - could not identify simulation: %s \n',simDisplayName);
end
fprintf('Setting stellar mass threshold to: %0.1e solar mass \n',massThresh);




%% select galaxies
M200c=double(fofs.Group_M_Crit200(subsInfo.hostFof+1)).*illUnits.massUnit;
R200c=double(fofs.Group_R_Crit200(subsInfo.hostFof+1));  %.*illUnits.lengthUnit;
T200c=double(fofs.Group_T_Crit200(subsInfo.hostFof+1));
K200c=mass_entropy_relation(M200c,'zred',illUnits.zred,'cosmo',cosmoStruct,...
    'delta',200,'rhotype','crit');
massAllGals=illustris.utils.get_stellar_mass(subs,'gal');
% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMaskBase=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'massTop',massThreshTop,'snap',snap,'gas','centrals');
galMask2=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'gas','centrals');

indxBase=find(galMaskBase);
indx2=find(galMask2);


nneib=illustris.utils.find_k_nearest_neighbor(subs.SubhaloPos(:,galMask2),1,'qp',subs.SubhaloPos(:,galMaskBase));
indx3=indx2(nneib.indx);
rnorm=max(R200c(indxBase),R200c(indx3));
isolatedMask=nneib.distance./rnorm>=10;

dwarfIndx=indxBase(isolatedMask);



%massHistStruct.mask=galMaskBase;
ids=dwarfIndx-1;
ids=shuffleArray(ids);
ids=ids(1:10);
massHistStruct.ids=ids;

%PropStruct.galMass=massAllGals;
%PropStruct.galMask=galMask;



%% generate values
%fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);

len2=length(ids);%sum(galMask);

% size of mass histogram
distLen=256;


%% initialize to zero
fprintf('Initializing output... \n');

% set component
compNames=["Gal", "CGMin", "CGMout", "CGM", "CGMoutskirt",  "CGMall", "Sub"];
paramNames=["Tcool",  "Temp", "Entropy", "Density"]; %"TcTff",
propNames=["MeanMW", "StdDevMW", "MassMedian", "MassQuantiles"];

% These values are not needed for the catalogs
% tcBin=[-1 0 1];
% tctffBin=[-1 0 1];
% tempBin=[5 6.5];
% entBin=(-2:1:2);
% densBin=(-4:1:-1);

% PropStruct.tcBin=tcBin;
% PropStruct.tctffBin=tctffBin;
% PropStruct.tempBin=tempBin;
% PropStruct.entBin=entBin;
% PropStruct.densBin=densBin;


%% initialize output

for fld=compNames  % different components of the SubHalo
    PropStruct.(fld).mask=int8(galMaskBase);
    for param=paramNames % physical parameters we look at
        for prop=propNames % type of metric produced
            pname=fld + param + prop;
            switch prop
                case 'MassQuantiles'
                    PropStruct.(fld).(pname)=zeros(4,len2);
                otherwise
                    PropStruct.(fld).(pname)=zeros(1,len2);
            end
        end
    end

    PropStruct.(fld).(fld+"GasMass")=zeros(1,len2);
    PropStruct.(fld).(fld+"SfrMass")=zeros(1,len2);
    PropStruct.(fld).(fld+"AvgDens")=zeros(1,len2);
end


qus=[0.1 0.25 0.5 0.75 0.9];
PropStruct.qants=qus([1 2 4 5]);


PropStruct.M200c=M200c(ids+1);
PropStruct.R200c=R200c(ids+1).*illUnits.lengthUnit;
PropStruct.T200c=T200c(ids+1);
PropStruct.K200c=K200c(ids+1);
PropStruct.ids=ids;
%% run over subind objects
fprintf(' *** Running over %s Galaxies *** \n',num2str(length(ids)));

cnt=0;
for id=ids

    structName="SubHalo_"+ num2str(id)

    % for following progress
    perCent=floor((cnt+1)/len2*100);

    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end

    % throw away unneeded objects
    %if galMaskBase(id+1)


    cnt=cnt+1;

    % load gas from in sub halo
    gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas');

    %indx=id+1;
    

    if gas.count==0
        continue
    end

    % identify Star-forming gas
    sfMask=gas.StarFormationRate>0;  %inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));

    % get density and temperature
    nDensity=double(gas.Density.*illUnits.numberDensityFactor); %in cm^-3
    gas=illustris.utils.addTemperature(gas);
    gas.Temperature(sfMask)=1000;
    %temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K

    %set mass
    mass=double(gas.Masses.*illUnits.massUnit); %in Solarmass

    % calculate cooling time
    if isfield(gas,'GFM_CoolingRate')
        tcool=double(illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate)); % cooling time in Gyr^-1
    else
        tcool=nan(size(mass));
    end

    % add entropy
    gas=illustris.utils.addEntropy( gas );


    % find distance from galaxy center
    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
    gasDist=sqrt( sum(double(gas.newCoord).^2,1)).*illUnits.lengthUnit;

    % find free-fall time
    % tff=recreate_TFF(gasDist,double(R200c(id+1)),...
    %     tffProfile.polyfit.a(id+1),...
    %     tffProfile.polyfit.b(id+1),...
    %     tffProfile.polyfit.c(id+1));
    % tctff=tcool./tff;

    % find important radii, all in physical units!
    rmax=max(gasDist);
    rhalfStar=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1)).*illUnits.lengthUnit; % stellar half mass radius
    rhalfGas=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1)).*illUnits.lengthUnit; % gas half mass radius
    rvir=double(R200c(id+1))*illUnits.lengthUnit;
    rgal=2.0.*rhalfStar;

    %PropStruct.rhalfStar=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1)); % stellar half mass radius
    tmpLim=log10([min(gas.Temperature(~sfMask)) max(gas.Temperature(~sfMask))]);
    dnsLim=log10([min(nDensity(~sfMask)) max(nDensity(~sfMask))]);
    %% go over components
    for fld=compNames

        switch fld
            case 'Gal'
                % for gas within the gal (2*r_half - stellar)
                distMask=gasDist<=rgal;

            case 'CGMin'
                % for gas within the cgm (2*r_half  - r_half gas)
                distMask=gasDist>rgal &  ...
                    gasDist<=rhalfGas;

            case 'CGMout'
                % gas from  half gas mass radius and rvir
                distMask=gasDist>rhalfGas  & ...
                    gasDist<=rvir;

            case 'CGMoutskirt'
                % gas from rvir outwards
                distMask=gasDist>rvir ;

            case 'CGM'
                % All CGM within R200,c
                distMask=gasDist>rgal & ...
                    gasDist<=rvir;

            case 'CGMall'
                % gas within 2* half stellar mass radius and edge
                distMask=gasDist>rgal ;

            case 'Sub'
                % all gas in the subfind object
                distMask=true(size(gasDist));
        end

        % build mask - non-sfr and within the ditance limitations
        mask=distMask & ~sfMask;

        if any(mask)
            % apply mask
            mm=mass(mask);
            tc=tcool(mask);
            %tcff=tctff(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            nDens=nDensity(mask);


            %% tcool and tc/tff
            % a small amount of cells have tc=0 set for positive cooling
            % rates (heating). We exclude them from the histogram and add
            % them to the last bin.
            tcMask=tc>0;
            if sum(tcMask)>0
                [mx, mhist, mus]=mk_mass_histogram(log10(tc(tcMask)),mm(tcMask),qus,distLen);
                %massDist(end)=massDist(end)+sum(mm(~tcMask));
                %[~,mxInd]=max(massDist);
                param="Tcool";
                %massHistStruct.(structName).(fld).
                PropStruct.(fld).(fld+param+"MeanMW")(cnt)=sum(mm(tcMask).*tc(tcMask))/sum(mm(tcMask));
                PropStruct.(fld).(fld+param+"StdDevMW")(cnt)=calc_standardDev(tc(tcMask),mm(tcMask));
                PropStruct.(fld).(fld+param+"MassMedian")(cnt)=10.^mus(3);
                PropStruct.(fld).(fld+param+"MassQuantiles")(:,cnt)=10.^mus([1 2 4 5]);

                massHistStruct.(structName).(fld).(param+"MassHist").hist=mhist;
                massHistStruct.(structName).(fld).(param+"MassHist").xAxis=mx;

                % [bird, binsize, xxlim,yylim]= histogram2d(log10(gas.Density),log10(gas.Temperature),...
                %     gas.Masses);



                % % tc / tff
                % [mx, mhist, mus]=mk_mass_histogram(log10(tcff(tcMask)),mm(tcMask),qus,distLen);
                % %massDist(end)=massDist(end)+sum(mm(~tcMask));
                % %[~,mxInd]=max(massDist);
                %
                % param='TcTff';
                % PropStruct.(fld).(strcat(fld,param,'MeanMW'))(indx)=sum(mm(tcMask).*tcff(tcMask))/sum(mm(tcMask));
                % PropStruct.(fld).(strcat(fld,param,'StdDevMW'))(indx)=calc_standardDev(tcff(tcMask),mm(tcMask));
                % PropStruct.(fld).(strcat(fld,param,'MassMedian'))(indx)=10.^mus(3);
                % PropStruct.(fld).(strcat(fld,param,'MassQuantiles'))(:,indx)=10.^mus([1 2 4 5]);
                %
            end

            %% temperature
            [mx, mhist, mus]=mk_mass_histogram(log10(tmp),mm,qus,distLen);
            %[~,mxInd]=max(massDist);
            param='Temp';
            PropStruct.(fld).(fld+param+"MeanMW")(cnt)=sum(mm.*tmp)/sum(mm);
            PropStruct.(fld).(fld+param+"StdDevMW")(cnt)=calc_standardDev(tmp,mm);
            PropStruct.(fld).(fld+param+"MassMedian")(cnt)=10.^mus(3);
            PropStruct.(fld).(fld+param+"MassQuantiles")(:,cnt)=10.^mus([1 2 4 5]);
            massHistStruct.(structName).(fld).(param+"MassHist").hist=mhist;
            massHistStruct.(structName).(fld).(param+"MassHist").xAxis=mx;

            %PropStruct.(fld).modeTemp(indx)=xx(mxInd);
            %
            %                 ll=length(tempBin)+1;
            %                 for i=1:ll
            %                     if i==1
            %                         msk=xx<=tempBin(i);
            %                     elseif i==ll
            %                         msk=xx>tempBin(i-1);
            %                     else
            %                         msk=xx>tempBin(i-1) & xx<=tempBin(i);
            %                     end
            %
            %                     PropStruct.(fld).massBinTemp(i,indx)=sum(massDist(msk));
            %                 end

            %% entropy
            [mx, mhist, mus]=mk_mass_histogram(log10(ent),mm,qus,distLen);
            %[~,mxInd]=max(massDist);
            param='Entropy';
            PropStruct.(fld).(fld+param+"MeanMW")(cnt)=sum(mm.*ent)/sum(mm);
            PropStruct.(fld).(fld+param+"StdDevMW")(cnt)=calc_standardDev(ent,mm);
            PropStruct.(fld).(fld+param+"MassMedian")(cnt)=10.^mus(3);
            PropStruct.(fld).(fld+param+"MassQuantiles")(:,cnt)=10.^mus([1 2 4 5]);
            massHistStruct.(structName).(fld).(param+"MassHist").hist=mhist;
            massHistStruct.(structName).(fld).(param+"MassHist").xAxis=mx;
            %PropStruct.(fld).modeEnt(indx)=xx(mxInd);

            %                 ll=length(entBin)+1;
            %                 for i=1:ll
            %                     if i==1
            %                         msk=xx<=entBin(i);
            %                     elseif i==ll
            %                         msk=xx>entBin(i-1);
            %                     else
            %                         msk=xx>entBin(i-1) & xx<=entBin(i);
            %                     end
            %
            %                     PropStruct.(fld).massBinEnt(i,indx)=sum(massDist(msk));
            %                 end

            %% number density
            [mx, mhist, mus]=mk_mass_histogram(log10(nDens),mm,qus,distLen);
            %[~,mxInd]=max(massDist);
            param='Density';
            PropStruct.(fld).(fld+param+"MeanMW")(cnt)=mean(nDens);
            PropStruct.(fld).(fld+param+"StdDevMW")(cnt)=std(nDens);
            PropStruct.(fld).(fld+param+"MassMedian")(cnt)=10.^mus(3);
            PropStruct.(fld).(fld+param+"MassQuantiles")(:,cnt)=10.^mus([1 2 4 5]);
            massHistStruct.(structName).(fld).(param+"MassHist").hist=mhist;
            massHistStruct.(structName).(fld).(param+"MassHist").xAxis=mx;
            %PropStruct.(fld).modeDensN(id+1)=xx(mxInd);

            %                 ll=length(densBin)+1;
            %                 for i=1:ll
            %                     if i==1
            %                         msk=xx<=densBin(i);
            %                     elseif i==ll
            %                         msk=xx>densBin(i-1);
            %                     else
            %                         msk=xx>densBin(i-1) & xx<=densBin(i);
            %                     end
            %
            %                     PropStruct.(fld).massBinDensN(i,id+1)=sum(massDist(msk));
            %                 end
            %


            [bird, binsize, xxlim,yylim]= histogram2d(log10(nDens),log10(tmp),...
                mm,'xlim',dnsLim,'ylim',tmpLim);
            massHistStruct.(structName).(fld).phaseDiagram.bird=bird;
            massHistStruct.(structName).(fld).phaseDiagram.binsize=binsize;
            massHistStruct.(structName).(fld).phaseDiagram.xlim=xxlim;
            massHistStruct.(structName).(fld).phaseDiagram.ylim=yylim;


            PropStruct.(fld).(fld+"GasMass")(cnt)=sum(mm); % only non-sf gas
            PropStruct.(fld).(fld+"SfrMass")(cnt)=sum(mass(distMask & sfMask));
            PropStruct.(fld).(fld+"AvgDens")(cnt)=sum(mm)/sum(mm./nDens); % in cm^-3
            PropStruct.(fld).([fld 'Sfr'])(cnt)=sum(gas.StarFormationRate(distMask & sfMask));
            %
            %                 PropStruct.(fld).cellNum(indx)=sum(mask);
            %
        end

    end

    %% build radial-parameter 2Dhistogrma
    % [bird, binsize, xxlim,yylim]= histogram2d(log10(nDensity(~sfMask)),log10(gas.Temperature(~sfMask)),...
    %     mass(~sfMask));
    % massHistStruct.(structName).phaseDiagram.bird=bird;
    % massHistStruct.(structName).phaseDiagram.binsize=binsize;
    % massHistStruct.(structName).phaseDiagram.xlim=xxlim;
    % massHistStruct.(structName).phaseDiagram.ylim=yylim;

    %rv=R200c(id+1);
    for p=paramNames
        switch p
            case "Tcool"
                par=tcool(~sfMask);
            case "Temp"
                par=gas.Temperature(~sfMask);
            case "Entropy"
                par=gas.Entropy(~sfMask);
            case "Density"
                par=nDensity(~sfMask);
        end

        [bird, binsize, xxlim,yylim]= histogram2d(log10(gasDist(~sfMask)./rvir),log10(par),...
            mass(~sfMask));
        massHistStruct.(structName).("rad" + p).bird=bird;
        massHistStruct.(structName).("rad" + p).binsize=binsize;
        massHistStruct.(structName).("rad" + p).xlim=xxlim;
        massHistStruct.(structName).("rad" + p).ylim=yylim;
    end
    massHistStruct.(structName).rgal=rgal;
    massHistStruct.(structName).rgas=rhalfGas;
    massHistStruct.(structName).R200c=rvir;
    %massHistStruct.(structName).R200c_simU=R200c;
    massHistStruct.(structName).rmax=rmax;
    massHistStruct.(structName).M200c=M200c(id+1);
    massHistStruct.(structName).T200c=T200c(id+1);
    massHistStruct.(structName).K200c=K200c(id+1);




    %end

end


%% save to catalog file
global DRACOFLAG
if DRACOFLAG



    % for fld=compNames
    % 
    %     outStruct=PropStruct.(fld);
    % 
    %     catName=sprintf('Subhalos_%s_PhysicalGasProperties',fld);
    % 
    %     folder=['PhysicalGasProperties/' simDisplayName];
    % 
    %     illustris.utils.write_catalog(outStruct,snap,'name',catName,...
    %         'path','default','folder','PhysicalGasProperties','v');
    % 
    % 
    %     %
    %     %         fname2=sprintf('gasProperties_massHistograms_snp%s_%s',num2str(snap),simDisplayName);
    %     %         save([DEFAULT_MATFILE_DIR '/' fname2],'massHist','-v7.3')
    %     %
    %     %         fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname2]);
    %     %
    % end

    fname=sprintf('gasProperties_isoDwarves_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'PropStruct','massHistStruct','-v7.3')

    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);

end
fprintf(' *** DONE!  *** \n');

%% save to postprocessed catalog

%global myPOSTPROCESSING
% PropStruct.inGal.mask=floor(galMask);
% PropStruct.inCGM.mask=floor(galMask);
% PropStruct.inSub.mask=floor(galMask);
% illustris.utils.write_catalog( PropStruct.inGal,99,'name','gasProps_inGal','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( PropStruct.inCGM,99,'name','gasProps_inCGM','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( PropStruct.inSub,99,'name','gasProps_inSub','folderName','gasProperties','verbose');
%



