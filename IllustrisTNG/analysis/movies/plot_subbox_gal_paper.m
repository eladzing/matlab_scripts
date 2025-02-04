%% perliminaries

%if ~skipFlag
    sim='100';
    
    bp=illustris.set_env(sim);
    %global BASEPATH;
    global DEFAULT_MATFILE_DIR
    global illUnits
    global cosmoStruct
    %global DEFAULT_MOVIE_DIR
    load([DEFAULT_MATFILE_DIR '/sub0_cgmList_TNG100.mat'])
    
    treeId=1;
    galID=tree(treeId).SubfindID(1); % 194946;
    
    %snaps=tree(treeId).SnapNum;
    %zredTNG=illustris.utils.snap2redshift(snaps);
    
    rLimitOuter=tree(treeId).SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,1); % gas half mass radius
    rLimitInner=tree(treeId).SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,1); % gas half mass radius
    
    bhMass=tree(treeId).SubhaloMassInRadType(illustris.partTypeNum('bh')+1,:).*illUnits.massUnit;
    stellarMass=tree(treeId).SubhaloMassInRadType(illustris.partTypeNum('star')+1,:).*illUnits.massUnit;
    m200c=tree(treeId).Group_M_Crit200.*illUnits.massUnit;
    ztree=illustris.utils.snap2redshift(tree(treeId).SnapNum);
    
    %% set subbox
    
    bpSub=illustris.set_env(sim,'sub0');
    
    %%read particles
    
    global subBox
    global simDisplayName
    snaps=0:subBox.Nsnaps-1;
    
    zred=illustris.utils.snap2redshift(snaps);
    times=redshift2time(zred);
    
    %     frameDir=[DEFAULT_MOVIE_DIR '/frames_id' num2str(galID) '_' simDisplayName];
    %     mkdir(frameDir)
    
    % identify centers
    load(sprintf('%s/cgmList_sub%i_posHistory_z4.01_%s.mat',DEFAULT_MATFILE_DIR,galID,simDisplayName));
    mask=subCenter(1,:)~=0;
    snp=subSnap(mask);
    tim=times.age(snp+1);
    %tim=fliplr(times.age(snp+1));
    
    cx=subCenter(1,mask);
    cy=subCenter(2,mask);
    cz=subCenter(3,mask);
    
    
    
    %% set times of interest
    zredStart=4;
    zredEnd=0;
    
    inds=find(zred<=zredStart & zred>=zredEnd);
    
    % interpolate centers
    cen(1,:)=interp1(tim,cx,times.age(inds));
    cen(2,:)=interp1(tim,cy,times.age(inds));
    cen(3,:)=interp1(tim,cz,times.age(inds));
    
    
    
    
    skipStep=10;
    m=true(size(inds));
    m(2:end-1)=mod(2:length(inds)-1,skipStep)==false;
    
%end


%% run over timesteps

fprintf('running over %i snaps \n',sum(m))


sparseInds=inds(m);
sparseCen=cen(:,m);
%cnt=0;
tic
for i= 1:length(sparseInds)
    
    
    
    indx=sparseInds(i);
    
    snap=snaps(indx);
    
    ztInd=find( zred(indx)<ztree,1,'first');
    
    mgal=mk_exponent_string(stellarMass(ztInd));
    mhalo=mk_exponent_string(m200c(ztInd));
    mBH=mk_exponent_string(bhMass(ztInd));
    
    msTag=['$M_\mathrm{gal}=' mgal '\,\mathrm{M_\odot}$'];
    mhTag=['$M_\mathrm{host}=' mhalo '\,\mathrm{M_\odot}$'];
    bhTag=['$M_\mathrm{BH}=' mBH '\,\mathrm{M_\odot}$'];
    
    fprintf('starting indx %i, snap %i zred %1.2f - ',i,snap,zred(indx));
    
    illustris.utils.set_illUnits(snap);
    
    
    stars=illustris.snapshot.loadSubset(bpSub,snap,illustris.partTypeNum('stars'),{'Coordinates','Masses','GFM_StellarFormationTime'});
    gas=illustris.snapshot.loadSubset(bpSub,snap,illustris.partTypeNum('gas'));
    gas=illustris.utils.addEntropy(gas);
    gas=illustris.utils.addTemperature(gas);
    gas=illustris.utils.addPressure(gas);
    
    center=sparseCen(:,i);
    
    for k=1:3
        stars.newCoord(k,:)=stars.Coordinates(k,:)-center(k);
        gas.newCoord(k,:)=gas.Coordinates(k,:)-center(k);
    end
    
    % remove unwanted particles
    
    %rrs=sqrt(sum(stars.newCoord.^2,1));
    %rrg=sqrt(sum(gas.newCoord.^2,1));
    
    rLim=1.05.*rLimitOuter;
    boxx=2*rLim;
    
    
    gasMask=abs(gas.newCoord(1,:))<=rLim & abs(gas.newCoord(2,:))<=rLim & abs(gas.newCoord(3,:))<=rLim ;
    starMask=abs(stars.newCoord(1,:))<=rLim & abs(stars.newCoord(2,:))<=rLim & abs(stars.newCoord(3,:))<=rLim ;
    starMask=starMask & stars.GFM_StellarFormationTime>=0; % remove wind particles
    %starMask=abs(stars.newCoord(1,:))<=rLimitOuter & abs(stars.newCoord(2,:))<=rLimitOuter & abs(stars.newCoord(3,:))<=rLimitOuter ;
    
    % find vcm
    
    vcm=sum(gas.Velocities(:,gasMask).*gas.Masses(gasMask),2)./sum(gas.Masses(gasMask));
    
    
    
    hf=figure;
    set(hf,'position',[32  68  2520  1260],'color','k','Visible','off');
    
    barProp.fontsize=22;
    %barProp.location='east';
    barProp.location='south';   %[0.785,0.135,0.028,0.765];
    barProp.color='w';
    barProp.tickLabelInterpreter='latex';
    barProp.axisLocation='in';
    
    ng=256;
    
    textM.str=msTag;%   ['$' mgalText '\,\mathrm{M_\odot}$'];
    textM.position=[-0.9*boxx/2 0.85*boxx/2 ];
    textM.color=[1 1 1];
    textM.fontsize=28;
    textM.Rotation=0;
    textH.str=mhTag;%   ['$' mgalText '\,\mathrm{M_\odot}$'];
    textH.position=[-0.9*boxx/2 0.85*boxx/2 ];
    textH.color=[0 0 0];
    textH.fontsize=28;
    textH.Rotation=0;
    textBH.str=bhTag;%   ['$' mgalText '\,\mathrm{M_\odot}$'];
    textBH.position=[-0.9*boxx/2 0.85*boxx/2 ];
    textBH.color=[0 0 0];
    textBH.fontsize=28;
    textBH.Rotation=0;
    
    
    %% stars
    newStarThresh=0.1;
    newStarTag='100\,\mathrm{Myr}';
    newStarColor=[0.09,0.63,0.86];
    
    scaleLength=100;
    scaleArw.start=[0.8*boxx/2-scaleLength 0.8*boxx/2];
    scaleArw.stop=[0.8*boxx/2 0.8*boxx/2];
    scaleArw.ends='none';
    scaleArw.color='w';
    
    
    textS.str=['$' num2str(scaleLength) '\,\mathrm{kpc}$'];
    textS.position=[0.8*boxx/2-scaleLength 0.86*boxx/2 ];
    textS.color=[1 1 1];
    textS.fontsize=22;
    textS.Rotation=0;
        
    textNS.str=['Young Stars $(<' newStarTag ')$'];
    textNS.position=[-0.9*boxx/2 0.74*boxx/2 ];
    textNS.color=newStarColor;
    textNS.fontsize=28;
    textNS.Rotation=0;
    
    text(1)=textS;
    text(2)=textM;
    text(3)=textNS;
    
    %%subplot(2,3,1)
    axStar=axes('Parent',hf,...
        'Position',[0.01 0.5  0.225 0.45],'color','w');
    %hold(axStar,'on');
    illustris.plots.mkmapStars('star',stars,'type','mass','ng',ng,'mask',starMask,'xz','box',boxx,'nn',3,...
        'figure',hf,'axes',axStar,'zoom',rLimitOuter,'clims',[4.5 9],'labels','no','arrow',scaleArw,'text',text,...
        'noticks','nogrid','barprop',barProp,'bartag',' ','newStars',newStarThresh,'newStarscolor',newStarColor);
    
    cmapStar=brewermap(256,'*Greys');
    %set(axStar,'XTickLabel','','YTickLabel','')
    ht=titlemine(' log Stellar Surface Density $[M_\odot/\mathrm{kpc^2}]$',26);set(ht,'color','w');
    
    %% temp
    %%subplot(2,3,2)
    barProp.color='k';
    axTemp=axes('Parent',hf,...
        'Position',[0.25 0.5  0.225 0.45]);
    illustris.plots.mkmapGas('gas',gas,'type','temp','xz','ng',ng,'mask',gasMask,'vfield','dilute',6,'vcm',vcm,'figure',hf,'axes',axTemp,...
        'box',boxx,'zoom',rLimitOuter, 'clims',[5 7.55],'streamdense',0,'labels','no','black',...
        'noticks','nogrid','barprop',barProp,'bartag',' ');
    
    set(axTemp,'XTickLabel','','YTickLabel','')
    cmapTemp=brewermap(256,'*Spectral');
    ht=titlemine('log Gas Temperature [K]',26);set(ht,'color','w');
     
    %% metallicity
    %subplot(2,3,5)
    axMet=axes('Parent',hf,...
        'Position',[0.75 0.5  0.225 0.45]);
    
    
    illustris.plots.mkmapGas('gas',gas,'type','ztot','xz','ng',ng,'mask',gasMask,'log','zoom',rLimitOuter,...
        'box',boxx,'clims',[-3.5 -1.45],'figure',hf,'axes',axMet,'labels','no','black',...
        'noticks','nogrid','barprop',barProp,'bartag',' ');
    set(axMet,'XTickLabel','','YTickLabel','')
    cmapMet=brewermap(256,'YlOrBr');
    
    ht=titlemine('log Gas Metallicity $[Z/Z_\odot]$',26);set(ht,'color','w');
    
    
     %% Pressure
    axPres=axes('Parent',hf,...
        'Position',[0.5 0.5  0.225 0.45]);
    
    text=textBH;
    illustris.plots.mkmapGas('gas',gas,'type','Pressure','xz','ng',ng,'mask',gasMask,'log','zoom',rLimitOuter,'text',text,...
        'box',boxx,'clims',[-4 -1.45],'figure',hf,'axes',axPres,'labels','no','black','noticks','nogrid','barprop',barProp,'bartag',' ');
    set(axPres,'XTickLabel','','YTickLabel','')
    cmapPres=brewermap(256,'PuRd');
    
    ht=titlemine('log Gas Pressure $[\mathrm{10^{-10}erg/cm^3}]$',26);
    %set(ht,'color','w');
    
    
    %% density
    
    if zred(indx)<1.5
        barProp.color='w';
        textH.color=[1 1 1];
    else
        barProp.color='k';
    end
    %barProp.location='South';
    
    %subplot(2,3,4)
    axDens=axes('Parent',hf,...
        'Position',[0.01 0.01 0.225 0.45]);
    text=textH;
    illustris.plots.mkmapGas('gas',gas,'type','ndensity','xz','ng',ng,'mask',gasMask,'figure',hf,'axes',axDens,...
        'box',boxx,'zoom',rLimitOuter, 'clims',[-4.5 -2],'labels','no','black','text',text,...
        'noticks','nogrid','barprop',barProp,'bartag',' ');
    
    %     cmapDens=brewermap(256,'OrRd');
    cmapDens=brewermap(256,'*YlGnBu');
    %cmapDens=brewermap(256,'YlOrRd');%YlGnBu');
    ht=titlemine('log Gas Density $[\mathrm{cm^{-3}}]$',26);set(ht,'color','w');
    
    %% entropy
    if zred(indx)<1.5
        barProp.color='w';
    else
        barProp.color='k';
    end
    %subplot(2,3,3)
    axEnt= axes('Parent',hf,...
        'Position',[0.25 0.01 0.225 0.45]);
    illustris.plots.mkmapGas('gas',gas,'type','entropy','xz','ng',ng,'mask',gasMask,'vfield','dilute',6,'vcm',vcm,'figure',hf,'axes',axEnt,...
        'box',boxx,'zoom',rLimitOuter,'clims',[0 3],'streamdense',0,'labels','no','black',...
        'noticks','nogrid','barprop',barProp,'bartag',' ');
    
    
    set(axEnt,'XTickLabel','','YTickLabel','')
    %cmapEnt=brewermap(256,'*BuPu');
    ht=titlemine('log Gas Entropy $[\mathrm{KeV\, cm^2}]$',26);set(ht,'color','w');
    %cmapEnt=brewermap(256,'*YlGnBu'); %*YlOrRd');
    cmapEnt=flipud(plasma(256));
    
  
    %% tcool
    %subplot(2,3,5)
    
           barProp.color='w';
    
    axtcool=axes('Parent',hf,...
        'Position',[0.5 0.01 0.225 0.45]);
    
    
    illustris.plots.mkmapGas('gas',gas,'type','tcool','xz','ng',ng,'mask',gasMask,'log','zoom',rLimitOuter,'clims',[-2 1.16],...
        'box',boxx,'figure',hf,'axes',axtcool,'labels','no','black','noticks','nogrid','barprop',barProp,'bartag',' ');
    set(axtcool,'XTickLabel','','YTickLabel','')
    cmapTc=flipud(viridis(256));
    ht=titlemine(' log Cooling Time $[\mathrm{Gyr}]$',26);set(ht,'color','w');
    
    
      %% vr
    %subplot(2,3,6)
    barProp.color='k';
    barProp.ticks=[-500 -250 0 250 500];
    axVr=axes('Parent',hf,...
        'Position',[0.75 0.01 0.225 0.45]);
    
    
    illustris.plots.mkmapGas('gas',gas,'type','vr','xz','ng',ng,'mask',gasMask,'vfield','dilute',6,'vcm',vcm,'figure',hf,'axes',axVr,...
        'box',boxx,'zoom',rLimitOuter,'clims',[-500 500],'streamdense',0,'labels','no','black',...
        'noticks','nogrid','barprop',barProp,'bartag',' ');
    % set(axVr,'XTickLabel','','YTickLabel','')
    cmapVr=brewermap(256,'*RdBu');
    ht=titlemine('Gas Radial Velocity $[\mathrm{km/sec}]$',26);set(ht,'color','w');
    
   
    %
    %     %% mach
    %     %subplot(2,3,5)
    %     axMach=axes('Parent',hf,...
    %         'Position',[0.5 0.01 0.225 0.45]);
    %
    %
    %     illustris.plots.mkmapGas('gas',gas,'type','mach','xz','ng',256,'mask',gasMask,'log','zoom',rLimitOuter,'clims',[0 1],...
    %         'box',boxx,'figure',hf,'axes',axMach,'labels','no','black')
    %     set(axMach,'XTickLabel','','YTickLabel','')
    %     cmapMach=brewermap(256,'*YlOrRd');
    %     cmapMach(1,:)=[0 0 0];
    %     ht=titlemine('Gas Mach Number',20);set(ht,'color','w');
    %
    
   
    
    
    %% set colors
    
    colormap(axStar,cmapStar);
    colormap(axTemp,cmapTemp);
    colormap(axDens,cmapDens);
    colormap(axEnt,cmapEnt);
    colormap(axVr,cmapVr);
    colormap(axtcool,cmapTc);
    colormap(axMet,cmapMet);
    colormap(axPres,cmapPres);
    
    age=redshift2time(zred(indx),'cosmo',cosmoStruct);
    
    zstr=sprintf('$z=%1.2f$, Age$= %1.2f$ Gyr',zred(indx),age.age);
    annotation(hf,'textbox',...
        [0.43 0.96 0.2 0.05],...
        'string',zstr,'Fontsize',30,...
        'Interpreter','latex',...
        'color',[1 1 1],'linestyle','none')
    
    annotation(hf,'textbox',...
        [0.93 0.95 0.2 0.05],...
        'string','E.\@ Zinger','Fontsize',26,...
        'Interpreter','latex',...
        'color',[1 1 1],'linestyle','none')
    
    
    clear barProp
    %    set(axStar,'colormap',cmapStar);
    %    set(axTemp,'colormap',cmapTemp);
    %    set(axDens,'colormap',cmapDens);
    %    set(axEnt,'colormap',cmapEnt);
    %    set(axVr,'colormap',cmapVr);
    %
    
    
    %illustris.plots.mkmapGas('gas',gas,'type','entropy','ng',256,'mask',gasMask,'box',0.5.*rLimitOuter,'avg')
    
    %cnt=cnt+1;
    
    F(i)=getframe(hf);
    %     zstr2=sprintf('%1.2f',zred(indx));
    %     fname=['snp' num2str(snap) '_z' zstr2];
    %     printout_fig(hf,fname,'dir',frameDir,'png','nopdf');
    %
    
%    printout_fig(hf,'testFig');
    close(hf)
    
    
    fprintf('completeted %s %% \n',num2str(100*i/length(sparseInds)))
 %  pause 
end
toc



mname=['/subMov_all_id' num2str(galID) '_z4_z0_' simDisplayName];
%mname2=[DEFAULT_MATFILE_DIR '/subMov_stars_id' num2str(galID) '_' simDisplayName];
save([DEFAULT_MATFILE_DIR '/' mname ],'F','-v7.3')
%T_MATFILE_DIR '/' mname ],'F','F2','-v7.3')
fprintf('done: %s',mname)

%galaxyMovieAssembly(F,mname,10);
%galaxyMovieAssembly(F2,mname);

