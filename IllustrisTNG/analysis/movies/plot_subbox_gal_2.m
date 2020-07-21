%% perliminaries

%if ~skipFlag
sim='100';

bp=illustris.set_env(sim);
%global BASEPATH;
global DEFAULT_MATFILE_DIR
global illUnits
load([DEFAULT_MATFILE_DIR '/sub0_cgmList_TNG100.mat'])

treeId=1;
galID=tree(treeId).SubfindID(1); % 194946;

%snaps=tree(treeId).SnapNum;
%zredTNG=illustris.utils.snap2redshift(snaps);

rLimitOuter=tree(treeId).SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,1); % gas half mass radius
rLimitInner=tree(treeId).SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,1); % gas half mass radius


%% set subbox

bpSub=illustris.set_env(sim,'sub0');

%%read particles

global subBox
global simDisplayName
snaps=0:subBox.Nsnaps-1;

zred=illustris.utils.snap2redshift(snaps);
times=redshift2time(zred);


% identify centers
load(sprintf('%s/cgmList_sub%i_posHistory_%s.mat',DEFAULT_MATFILE_DIR,galID,simDisplayName));
mask=subCenter(1,:)~=0;
snp=subSnap(mask);
tim=times.age(snp+1);
%tim=fliplr(times.age(snp+1));

cx=subCenter(1,mask);
cy=subCenter(2,mask);
cz=subCenter(3,mask);



%% set times of interest
zredStart=2;
zredEnd=0.5;

inds=find(zred<=zredStart & zred>=zredEnd);

% interpolate centers
cen(1,:)=interp1(tim,cx,times.age(inds));
cen(2,:)=interp1(tim,cy,times.age(inds));
cen(3,:)=interp1(tim,cz,times.age(inds));





skipStep=5;
m=true(size(inds));
m(2:end-1)=mod(2:length(inds)-1,skipStep)==0;
%end


%% run over timesteps




sparseInds=inds(m);
sparseCen=cen(:,m);

fprintf('running over %i snaps from z=%s to z=%s \n',sum(m),num2str(zred(sparseInds(1))),num2str(zred(sparseInds(end))));
%cnt=0;
tic
for i=1:length(sparseInds)
    
    
    
    indx=sparseInds(i);
    
    snap=snaps(indx);
    
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
    set(hf,'position',[32  68  2520  1260],'Visible','off');
    
    %% stars
    scaleArw.start=[-300 -300];
    scaleArw.stop=[-200 -300];
    scaleArw.width=2;
    scaleArw.color='w';
    scaleArw.ends='none';
    
    text.str='$100\,\mathrm{ckpc}$';
    text.pos=[-320 -240 200 20];
    text.color=[1 1 1];
    
    %%subplot(2,3,1)
    axStar=axes('Parent',hf,...
        'Position',[0.01 0.48  0.23 0.46]);
    %hold(axStar,'on');
    illustris.plots.mkmapStars('star',stars,'type','mass','ng',256,'mask',starMask,'xz','box',boxx,...
        'figure',hf,'axes',axStar,'zoom',rLimitOuter,'clims',[5 9],'labels','no','arrow',scaleArw,'text',text);
    cmapStar=brewermap(256,'*Greys');
    set(axStar,'XTickLabel','','YTickLabel','')
    titlemine('Stellar Surface Density')
    
    %% temp
    %%subplot(2,3,2)
    axTemp=axes('Parent',hf,...
        'Position',[0.25 0.48  0.23 0.46]);
    
    illustris.plots.mkmapGas('gas',gas,'type','temp','xz','ng',256,'mask',gasMask,'vfield','dilute',6,'vcm',vcm,'figure',hf,'axes',axTemp,...
        'box',boxx,'zoom',rLimitOuter, 'clims',[5 7.5],'streamdense',0,'labels','no')
    set(axTemp,'XTickLabel','','YTickLabel','')
    cmapTemp=brewermap(256,'*Spectral');
    titlemine('Gas Temperature')
    
    %% density
    %subplot(2,3,4)
    axDens=axes('Parent',hf,...
        'Position',[0.01 0.01 0.23 0.46]);
    
    illustris.plots.mkmapGas('gas',gas,'type','ndensity','xz','ng',256,'mask',gasMask,'figure',hf,'axes',axDens,...
        'box',boxx,'zoom',rLimitOuter, 'clims',[-4.5 -2],'labels','no')
    set(axDens,'XTickLabel','','YTickLabel','')
    %     cmapDens=brewermap(256,'OrRd');
    cmapDens=brewermap(256,'YlOrRd');%YlGnBu');
    titlemine('Gas Density')
    
    %% entropy
    %subplot(2,3,3)
    axEnt= axes('Parent',hf,...
        'Position',[0.25 0.01 0.23 0.46]);
    illustris.plots.mkmapGas('gas',gas,'type','entropy','xz','ng',256,'mask',gasMask,'vfield','dilute',6,'vcm',vcm,'figure',hf,'axes',axEnt,...
        'box',boxx,'zoom',rLimitOuter,'clims',[1 3],'streamdense',0,'labels','no');
    set(axEnt,'XTickLabel','','YTickLabel','')
    %cmapEnt=brewermap(256,'*BuPu');
    titlemine('Gas Entropy')
    cmapEnt=brewermap(256,'*YlGnBu'); %*YlOrRd');
    
    
    %% vr
    %subplot(2,3,6)
    axVr=axes('Parent',hf,...
        'Position',[0.75 0.01 0.23 0.46]);
    
    
    illustris.plots.mkmapGas('gas',gas,'type','vr','xz','ng',256,'mask',gasMask,'vfield','dilute',6,'vcm',vcm,'figure',hf,'axes',axVr,...
        'box',boxx,'zoom',rLimitOuter,'clims',[-500 500],'streamdense',0,'labels','no')
    set(axVr,'XTickLabel','','YTickLabel','')
    cmapVr=brewermap(256,'*RdBu');
    titlemine('Gas Radial Velocity')
    
    %% metallicity
    %subplot(2,3,5)
    axMet=axes('Parent',hf,...
        'Position',[0.75 0.48  0.23 0.46]);
    
    
    illustris.plots.mkmapGas('gas',gas,'type','ztot','xz','ng',256,'mask',gasMask,'log','zoom',rLimitOuter,...
        'box',boxx,'clims',[-3.5 -1.5],'figure',hf,'axes',axMet,'labels','no')
    set(axMet,'XTickLabel','','YTickLabel','')
    cmapMet=brewermap(256,'YlOrBr');
    
    titlemine('Gas Metallicity')
    
    %% mach
    %subplot(2,3,5)
    axMach=axes('Parent',hf,...
        'Position',[0.5 0.01 0.23 0.46]);
    
    
    illustris.plots.mkmapGas('gas',gas,'type','mach','xz','ng',256,'mask',gasMask,'log','zoom',rLimitOuter,'clims',[0 1],...
        'box',boxx,'figure',hf,'axes',axMach,'labels','no')
    set(axMach,'XTickLabel','','YTickLabel','')
    cmapMach=brewermap(256,'*YlOrRd');
    cmapMach(1,:)=[0 0 0];
    titlemine('Gas Mach Number')
    
    
    %% Pressure
    axPres=axes('Parent',hf,...
        'Position',[0.5 0.48  0.23 0.46]);
    
    
    illustris.plots.mkmapGas('gas',gas,'type','Pressure','xz','ng',256,'mask',gasMask,'log','zoom',rLimitOuter,...
        'box',boxx,'clims',[-4 -1.5],'figure',hf,'axes',axPres,'labels','no')
    set(axPres,'XTickLabel','','YTickLabel','')
    cmapPres=brewermap(256,'PuRd');
    
    titlemine('Gas Pressure')
    
    
    %% set colors
    
    colormap(axStar,cmapStar);
    colormap(axTemp,cmapTemp);
    colormap(axDens,cmapDens);
    colormap(axEnt,cmapEnt);
    colormap(axVr,cmapVr);
    colormap(axMach,cmapMach);
    colormap(axMet,cmapMet);
    colormap(axPres,cmapPres);
    
    zstr=sprintf('$z=%1.2f $',zred(indx));
    annotation(hf,'textbox',...
        [0.46 0.95 0.08 0.05],...
        'string',zstr,'Fontsize',25,...
        'Interpreter','latex',...
        'color',[0 0 0],'linestyle','none')
    
    
    
    %    set(axStar,'colormap',cmapStar);
    %    set(axTemp,'colormap',cmapTemp);
    
     
    %    set(axEnt,'colormap',cmapEnt);
    %    set(axVr,'colormap',cmapVr);
    %
    
    
    %illustris.plots.mkmapGas('gas',gas,'type','entropy','ng',256,'mask',gasMask,'box',0.5.*rLimitOuter,'avg')
    
    %cnt=cnt+1;
    F(i)=getframe(hf);
    close(hf)
    
    
    fprintf('completeted %s %% \n',num2str(100*i/length(sparseInds)))
    
end
toc



mname=['/subMov_id' num2str(galID) '_z2-0.5_' simDisplayName '.mat'];
%mname2=[DEFAULT_MATFILE_DIR '/subMov_stars_id' num2str(galID) '_' simDisplayName];
save([DEFAULT_MATFILE_DIR '/' mname ],'F','-v7.3')
%T_MATFILE_DIR '/' mname ],'F','F2','-v7.3')
fprintf('done')

%galaxyMovieAssembly(F,mname,10);
%galaxyMovieAssembly(F2,mname);

