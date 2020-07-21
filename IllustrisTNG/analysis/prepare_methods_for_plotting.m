%% plot method exploration compared to tng100

method_list={'0000','2201','3000','3101','3102','3103','3104',...
    '3301','3302','3403','3404','3801','3802','2101','2302','0030'};



%methodDir='/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/methods/';
yl=[-3 3];
xl=[9 12.5];
snap35=4;

for i=1:length(method_list)
    
    methodNum=method_list{i};
    bp=illustris.set_env(35,methodNum);
    
    if i==1
        global DEFAULT_MATFILE_DIR
        
        methodDir=[DEFAULT_MATFILE_DIR '/methods'];
    end
    
    
    % load fofs and subs
    fofs35=illustris.groupcat.loadHalos(bp,snap35);
    subs35=illustris.groupcat.loadSubhalos(bp,snap35);
    
    fofs35=illustris.utils.addTvirFofs(fofs35);
    
    
    % load gas properties
    
    name=[methodDir '/gasProperties_snp4_TNG35_' methodNum '.mat'];
    load(name);
    
    subsInfo35 = illustris.infrastructure.build_sub_fof_connection(subs35,fofs35);
    centralMask35= subsInfo35.isCentral(tCoolStruct.galMask);
    
    tvir=fofs35.Group_T_Mean200(subsInfo35.hostFof+1);
    tvir=tvir(tCoolStruct.galMask);
    
    
    galMass=tCoolStruct.galMass(tCoolStruct.galMask);  % galaxy stellar mass
    %sfr=subs35.SubhaloSFRinRad(tCoolStruct.galMask);  % sfr in galaxy
    if ~strcmp(methodNum,'0030')
        ssfr=illustris.utils.calc_ssfr(subs35);% sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));
        ssfr=ssfr(tCoolStruct.galMask);
    else
        ssfr=zeros(size(tCoolStruct.galMask));
    end

% ssfrThresh=1e-11;
% qMask=ssfr<=ssfrThresh;


    
    
    %% gal 
    
    galTemp=tCoolStruct.inGal.meanTempMW(1,tCoolStruct.galMask)./tvir;
    galEnt=tCoolStruct.inGal.meanEntMW(1,tCoolStruct.galMask);
    galTc=tCoolStruct.inGal.meanTcMW(1,tCoolStruct.galMask);
    galDens=tCoolStruct.inGal.meanDensN(1,tCoolStruct.galMask);
    fgs=(tCoolStruct.inGal.gasMass(tCoolStruct.galMask)+...
        tCoolStruct.inGal.sfrMass(tCoolStruct.galMask))./galMass;
    fg=fgs./(1+fgs);
    
    galaxyMask=centralMask35 & galEnt>0;
    
    
    methodStruct.(['meth' methodNum]).galMass=log10(galMass(galaxyMask));
    methodStruct.(['meth' methodNum]).galSsfr=log10(ssfr(galaxyMask));
    methodStruct.(['meth' methodNum]).fgs=log10(fgs(galaxyMask));
    
    methodStruct.(['meth' methodNum]).galTemp=log10(galTemp(galaxyMask));
    methodStruct.(['meth' methodNum]).galEnt=log10(galEnt(galaxyMask));
    methodStruct.(['meth' methodNum]).galTc=log10(galTc(galaxyMask));
    methodStruct.(['meth' methodNum]).galDens=log10(galDens(galaxyMask));
       
    %% cgm 
       
    cgmTemp=tCoolStruct.inCGM.meanTempMW(1,tCoolStruct.galMask)./tvir;
    cgmEnt=tCoolStruct.inCGM.meanEntMW(1,tCoolStruct.galMask);
    cgmTc=tCoolStruct.inCGM.meanTcMW(1,tCoolStruct.galMask);
    cgmDens=tCoolStruct.inCGM.meanDensN(1,tCoolStruct.galMask);
    
    galaxyMask=centralMask35 & cgmEnt>0;
    
    
    methodStruct.(['meth' methodNum]).cgmMass=log10(galMass(galaxyMask));
    methodStruct.(['meth' methodNum]).cgmTemp=log10(cgmTemp(galaxyMask));
    methodStruct.(['meth' methodNum]).cgmEnt=log10(cgmEnt(galaxyMask));
    methodStruct.(['meth' methodNum]).cgmTc=log10(cgmTc(galaxyMask));
    methodStruct.(['meth' methodNum]).cgmDens=log10(cgmDens(galaxyMask));
    methodStruct.(['meth' methodNum]).cgmSsfr=log10(ssfr(galaxyMask));
    
    %% out 
    
    outTemp=tCoolStruct.inOut.meanTempMW(1,tCoolStruct.galMask)./tvir;
    outEnt=tCoolStruct.inOut.meanEntMW(1,tCoolStruct.galMask);
    outTc=tCoolStruct.inOut.meanTcMW(1,tCoolStruct.galMask);
    outDens=tCoolStruct.inOut.meanDensN(1,tCoolStruct.galMask);
        
    galaxyMask=centralMask35 & outEnt>0;
    
    methodStruct.(['meth' methodNum]).outMass=log10(galMass(galaxyMask));
    methodStruct.(['meth' methodNum]).outTemp=log10(outTemp(galaxyMask));
    methodStruct.(['meth' methodNum]).outEnt=log10(outEnt(galaxyMask));
    methodStruct.(['meth' methodNum]).outTc=log10(outTc(galaxyMask));
    methodStruct.(['meth' methodNum]).outDens=log10(outDens(galaxyMask));
    methodStruct.(['meth' methodNum]).outSsfr=log10(ssfr(galaxyMask));
    
    
end


fname=[DEFAULT_MATFILE_DIR '/methodStructure_snp4.mat'];
save(fname,'methodStruct','-v7.3')
    
