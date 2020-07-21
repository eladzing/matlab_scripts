global illUnits
global BASEPATH
%global fullSnapMask
global cosmoStruct

%rv=tree.Group_R_Crit200(1:15);
fofID=tree.SubhaloGrNr(1:15);


for i=1:15
    fprintf('i= %i \n',i);
    snap=tree.SnapNum(i);
    
    illustris.utils.set_illUnits(snap) % set the simulation units for the snapshot.
    
    rv(i)=tree.Group_R_Crit200(i).*illUnits.lengthUnit;
    [~, ~, tv(i), ~]=calculate_virials('mvir',tree.Group_M_Crit200(i)*illUnits.massUnit,'cosmo',cosmoStruct,'delv',200,'crit');
    rh(i)=tree.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,i).*illUnits.lengthUnit;

    gasFields={'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses','Velocities','StarFormationRate'};
    
    
    gasH=illustris.snapshot.loadHalo(BASEPATH,snap,fofID(i), 'gas',gasFields);
    mask=gasH.StarFormationRate==0;
    
    %nDensity=gasH.Density.*illUnits.numberDensityFactor; %in cm^-3
    mass=gasH.Masses.*illUnits.massUnit; %in Solarmreass
    
    
    
    gasH=illustris.utils.addTemperature(gasH);
    %gasH=illustris.utils.addCoolingTime(gasH);
    gasH=illustris.utils.addEntropy(gasH);
    gasH=illustris.utils.addPressure(gasH);
    
    gasH.newCoord = illustris.utils.centerObject(gasH.Coordinates,tree.GroupPos(:,i)).*illUnits.lengthUnit;
    
    
    gasH=illustris.utils.addCellRadius(gasH);
    
    %% Temperature
    resT(i) = mk_radial_profile_cells(gasH.newCoord,gasH.CellRadius,gasH.Temperature,'wt',mass,'type','intensive','nb',300,'mask',mask);
    
    %% Entropy
    
    %% Pressure
    resE(i) = mk_radial_profile_cells(gasH.newCoord,gasH.CellRadius,gasH.Entropy,'wt',mass,'type','intensive','nb',300,'mask',mask);
    %% Density
    resP1(i) = mk_radial_profile_cells(gasH.newCoord,gasH.CellRadius,gasH.Pressure,'wt',mass,'type','intensive','nb',300,'mask',mask);
    resP2(i) = mk_radial_profile_cells(gasH.newCoord,gasH.CellRadius,gasH.Pressure,'wt',ones(size(mass)),'type','intensive','nb',300,'mask',mask);
    resP3(i) = mk_radial_profile_cells(gasH.newCoord,gasH.CellRadius,gasH.Pressure,'wt',gasH.CellRadius,'type','intensive','nb',300,'mask',mask);
    figrue
    %% Cooling Time
    %restc(i) = mk_radial_profile_cells(gasH.newCoord,gasH.CellRadius,gasH.CoolingTime,'wt',mass,'type','intensive','log','nb',300);
    
end

