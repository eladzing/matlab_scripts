%% plot images for galaxies
simName='100';
snap=99;
bp=illustris.set_env(simName);

%global simDisplayName
global illUnits
global DEFAULT_PRINTOUT_DIR

printoutDir=[DEFAULT_PRINTOUT_DIR '/maps'];
close all
if perlimFlag
    %load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/fofs_subs_TNG100_z0.mat')
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
    
    
    %     massThresh=10^9; % threshold for *stellar* mass
    
    %     ssfr=illustris.utils.calc_ssfr(subs);
    
    %     global illUnits
    %
    %     massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
    
    %     galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals','hasGas');
    
end

cmapStar=brewermap(256,'*Greys');
cmapTemp=brewermap(256,'*Spectral');
cmapDens=brewermap(256,'*YlGnBu');
cmapEnt=brewermap(256,'PuRd');
cmapTc=brewermap(256,'YlOrRd');
cmapPress=brewermap(256,'YlOrBr');

caxTemp=[4.5 6.5];
caxEnt=[-0.5 2.5];
caxDens=[-5.5 -2];
caxTc=[-1 1.5];
caxPress=[-6 -3];

% 511214 460526 455730 488174
idList=[ 488841 490053 497088 499161 534514 ...
    566543 552112 570763 587022 595189 637717 350183 390859 405568 426575 ...
    438038 447153 452651 462077 472174 501526 511005 588650 594081 616552 ...
    627692 653913 658143 678411];


for ii=1:length(idList)
    idGal=idList(ii);
    proj='all';
    
    idHost=subsInfo.hostFof(idGal+1);
    
    printTag=['id' num2str(idGal) '_cgmFoc'];
    
    
    %baseName=['galMap_id%i_' proj '%s_snp' num2str(snap) '-' simDisplayName];
    
    %
    %         idGal=indx(i)-1;
    %
    galMass=double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,idGal+1).*illUnits.massUnit); % stellar mass within 2*rhalf
    %
    %massAllGals(idGal+1);
    
    %
    
    hostMass=fofs.Group_M_Crit200(idHost+1).*illUnits.massUnit;
    %
    
    center=subs.SubhaloPos(:,idGal+1).*illUnits.lengthUnit;
    
    rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,idGal+1).*illUnits.lengthUnit;
    rhalfStars=subs.SubhaloHalfmassRadType(illustris.partTypeNum('star')+1,idGal+1).*illUnits.lengthUnit;
    r200=fofs.Group_R_Crit200(idHost+1).*illUnits.lengthUnit;
    r500=fofs.Group_R_Crit500(idHost+1).*illUnits.lengthUnit;
    
    rLimitOuter=1.5*rhalfGas;
    rLim=1.05.*rLimitOuter;
    boxx=2*rLim;
    
    circ(1).radius=2.*rhalfStars;
    circ(1).width=1.5;
    circ(1).color=[0 0 0];
    circ(2).radius=rhalfGas;
    circ(2).width=1.5;
    circ(2).color=[0 0 0];
    circ(3).radius=r200;
    circ(3).width=1.5;
    circ(3).color=[1 1 1];
    circ(3).type=':';
    % circ(4).radius=r200;
    % circ(4).width=1.5;
    
    
    
    %% load  and  plot star data
    if readFlag
        starFields={'Coordinates','Masses','GFM_StellarFormationTime'};
        stars=illustris.snapshot.loadSubhalo(bp,snap,idGal,'stars',starFields);
        stars.newCoord=illustris.utils.centerObject(stars.Coordinates.*illUnits.lengthUnit,center);
        starMask=true(size(stars.newCoord(1,:)));    %   abs(stars.newCoord(1,:))<=rLim & abs(stars.newCoord(2,:))<=rLim & abs(stars.newCoord(3,:))<=rLim ;
        starMask=starMask & stars.GFM_StellarFormationTime>=0; % remove wind particles
        
        
        %'figure',hf,'axes',axStar,'zoom',rLimitOuter,'clims',[5 9],'labels','no','arrow',scaleArw,'text',text);
        %
        
        % load gas data
        gasFields={'Coordinates','Masses','Density','ElectronAbundance'...
            'EnergyDissipation','GFM_CoolingRate','GFM_Metallicity',...
            'InternalEnergy','Machnumber','MagneticField','NeutralHydrogenAbundance',...
            'Potential','StarFormationRate','Velocities'};
        
        gas=illustris.snapshot.loadHalo(bp,snap,idHost,'gas',gasFields);
        %gas=illustris.snapshot.loadSubhalo(bp,snap,idGal,'gas',gasFields);
        
        
        gas.newCoord=illustris.utils.centerObject(gas.Coordinates.*illUnits.lengthUnit,center);
        
        gas=illustris.utils.addEntropy(gas);
        gas=illustris.utils.addTemperature(gas);
        gas=illustris.utils.addPressure(gas);
        
        gasMask=true(size(gas.newCoord(1,:)));
        %gasMask=abs(gas.newCoord(1,:))<=rLim & abs(gas.newCoord(2,:))<=rLim & abs(gas.newCoord(3,:))<=rLim ;
        
    end
    %vcm=sum(gas.Velocities(:,gasMask).*gas.Masses(gasMask),2)./sum(gas.Masses(gasMask));
    
    %% plot
    ng=256;
    %titleTag=['gal ID: ' num2str(idGal)];
    
    mgalText=mk_exponent_string(galMass);
    mhostText=mk_exponent_string(hostMass);
    
    titleTag=['$' mgalText ',\,'  mhostText '\,\mathrm{M_\odot}$'];
    
    % ssfrString=mk_exponent_string(ssfr(idGal+1));
    % titleTag=['$M_\mathrm{s}=' mgalText' ',\, M_\mathrm{h}=' mhostText(2:end) '\,[\mathrm{M_\odot}]$']
    % if ssfr(idGal+1)>=1e-11
    %     qTag='SF';
    % else
    %     qTag='Q';
    % end
    
    
    
    
    thickness=1*rhalfGas;
    
    % run over projection
    
    
    %             hf=figure;
    %             set(hf,'position',[ 25          18        1548         967],'color','k');
    %
    
    %% stars
    
    %             axStar=axes('Parent',hf,...
    %                 'Position',[0.09 0.53  0.26 0.35],'color','w');
    %
    
    illustris.plots.mkmapStars('star',stars,'type','mass','ng',ng,proj,...
        'mask',starMask,'box',boxx,'cmap',cmapStar,...
        'zoom',rLimitOuter,'clims',[5 9],'circ',circ,'print',printTag,'savefig','printoutdir',printoutDir);
    
    %ht=titlemine('Stars',16);set(ht,'color','w');
    %             fname=sprintf(baseName,idGal,qTag,projTag{j},snap);
    %             printout_fig(hf,fname,'subdir','mapSurvey','nopdf')
    %
    %% density
    %subplot(2,3,4)
    %             axDens=axes('Parent',hf,...
    %                 'Position',[0.09 0.1 0.26 0.35]);
    illustris.plots.mkmapGas('gas',gas,'type','ndensity','ng',ng,proj,...
        'mask',gasMask,'circ',circ,'clims',caxDens,'cmap',cmapDens,...%'vfield','dilute',6,'vcm',vcm,...
        'box',boxx,'zoom',rLimitOuter,'thick',thickness,...
        'black','title',titleTag,'print',printTag,'savefig','printoutdir',printoutDir);  %,...%'streamdense',0,...
    
    %ht=titlemine('Density',16);set(ht,'color','w');
    
    
    %% entropy
    %             axEnt=axes('Parent',hf,...
    %                 'Position',[0.4 0.53  0.26 0.35]);
    %
    illustris.plots.mkmapGas('gas',gas,'type','entropy','ng',ng,proj,...
        'mask',gasMask,'circ',circ,'clims',caxEnt,'cmap',cmapEnt,...
        'box',boxx,'zoom',rLimitOuter,'thick',thickness,...
        'black','title',titleTag,'print',printTag,'savefig','printoutdir',printoutDir);
    
    %ht=titlemine('Entropy',16);set(ht,'color','w');
    
    %% temp
    %%subplot(2,3,2)
    %             axTemp=axes('Parent',hf,...
    %                 'Position',[0.4 0.1  0.26 0.35]);
    
    illustris.plots.mkmapGas('gas',gas,'type','temp','ng',ng,proj,...
        'mask',gasMask,'circ',circ,'clims',caxTemp,'cmap',cmapTemp,...%'vfield','dilute',6,'vcm',vcm,...
        'box',boxx,'zoom',rLimitOuter,'thick',thickness,...
        'black','title',titleTag,'print',printTag,'savefig','printoutdir',printoutDir); %,...%'streamdense',0,...
    
    %ht=titlemine('Temperature',16);set(ht,'color','w');
    
    
    %% Pressure
    %             axPres=axes('Parent',hf,...
    %                 'Position',[0.7 0.53  0.26 0.35]);
    
    illustris.plots.mkmapGas('gas',gas,'type','pressure','ng',ng,proj,...
        'mask',gasMask,'circ',circ,'cmap',cmapPress,'clims',caxPress,...%'vfield','dilute',6,'vcm',vcm,...
        'box',boxx,'zoom',rLimitOuter,'thick',thickness,...
        'black','title',titleTag,'print',printTag,'savefig','printoutdir',printoutDir); %,...%'streamdense',0,...
    
    %ht=titlemine('Pressure',16);set(ht,'color','w');
    
    %% tcool
    %             axTcool=axes('Parent',hf,...
    %                 'Position',[0.7 0.1 0.26 0.35]);
    illustris.plots.mkmapGas('gas',gas,'type','tcool','ng',ng,proj,...
        'mask',gasMask,'circ',circ,'clims',caxTc,'cmap',cmapTc,...%'vfield','dilute',6,'vcm',vcm,...
        'box',boxx,'zoom',rLimitOuter,'thick',thickness,'nobackground',...
        'black','title',titleTag,'print',printTag,'savefig','printoutdir',printoutDir); %,...%'streamdense',0,...
    
    fprintf('completed %i of %i \n',ii,length(idList))
   
end
%ht=titlemine('Cooling Time',16);set(ht,'color','w');



%% set colors

%             colormap(axStar,cmapStar);
%             colormap(axTemp,cmapTemp);
%             colormap(axDens,cmapDens);
%             colormap(axEnt,cmapEnt);
%             colormap(axTcool,cmapTc);
%             colormap(axPres,cmapPress);

%
%             zstr=titleTag;  %sprintf('$z=%1.2f $',zred(indx));
%             annotation(hf,'textbox',...
%                 [0.18 0.95 0.8 0.05],...
%                 'string',zstr,'Fontsize',20,...
%                 'Interpreter','latex',...
%                 'color',[1 1 1],'linestyle','none')




