%% perliminaries

if ~skipFlag
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
    zredStart=1;
    zredEnd=0;
    
    inds=find(zred<=zredStart & zred>=zredEnd);
    
    % interpolate centers
    cen(1,:)=interp1(tim,cx,times.age(inds));
    cen(2,:)=interp1(tim,cy,times.age(inds));
    cen(3,:)=interp1(tim,cz,times.age(inds));
    
    
    
    
    
    skipStep=20;
    m=true(size(inds));
    m(2:end-1)=mod(2:length(inds)-1,skipStep)==0;
end


%% run over timesteps


sparseInds=inds(m);
sparseCen=cen(:,m);
%cnt=0;
for i=1:length(sparseInds)
    tic
    
    
    indx=sparseInds(i);
    
    snap=snaps(indx);
    illustris.utils.set_illUnits(snap);
    
    
    stars=illustris.snapshot.loadSubset(bpSub,snap,illustris.partTypeNum('stars'),{'Coordinates','Masses'});
    gas=illustris.snapshot.loadSubset(bpSub,snap,illustris.partTypeNum('gas'));
    %gas=illustris.utils.addEntropy(gas);
    gas=illustris.utils.addTemperature(gas);
    %gas=illustris.utils.addPressure(gas);
    
    center=sparseCen(:,i);
    
    for k=1:3
        stars.newCoord(k,:)=stars.Coordinates(k,:)-center(k);
        gas.newCoord(k,:)=gas.Coordinates(k,:)-center(k);
    end
    
    % remove unwanted particles
    
    %rrs=sqrt(sum(stars.newCoord.^2,1));
    %rrg=sqrt(sum(gas.newCoord.^2,1));
    
    gasMask=abs(gas.newCoord(1,:))<=rLimitOuter & abs(gas.newCoord(2,:))<=rLimitOuter & abs(gas.newCoord(3,:))<=rLimitOuter ;
    starMask=abs(stars.newCoord(1,:))<=rLimitOuter & abs(stars.newCoord(2,:))<=rLimitOuter & abs(stars.newCoord(3,:))<=rLimitOuter ;
    %starMask=abs(stars.newCoord(1,:))<=rLimitOuter & abs(stars.newCoord(2,:))<=rLimitOuter & abs(stars.newCoord(3,:))<=rLimitOuter ;
    
    % find vcm
    
    vcm=sum(gas.Velocities(:,gasMask).*gas.Masses(gasMask),2)./sum(gas.Masses(gasMask));
    
    
    
     hf=figure;
    set(hf,'position',[32  68  2504  1261]);%76 
    hf=figure('Visible','off');
    
    illustris.plots.mkmapGas('gas',gas,'type','ndensity','xz','ng',256,'mask',gasMask,'vfield','dilute',8,'box',2.*rLimitOuter,'vcm',vcm,'figure',hf,...
        'clims',[-5 -1])
    
    hf2=figure('Visible','off');
    illustris.plots.mkmapGas('gas',gas,'type','ndensity','xz','ng',256,'mask',gasMask,'vfield','dilute',8,'box',2.*rLimitOuter,'vcm',vcm,'figure',hf,...
        'clims',[-5 -1])
    
    %illustris.plots.mkmapStars('star',stars,'type','mass','ng',256,'mask',starMask,'figure',hf2);
    
    
    %illustris.plots.mkmapGas('gas',gas,'type','entropy','ng',256,'mask',gasMask,'box',0.5.*rLimitOuter,'avg')
    
    %cnt=cnt+1;
    F(i)=getframe(hf);
    close(hf)
    
    F2(i)=getframe(hf2);
    close(hf2)
    
    fprintf('completeted %s %% ',num2str(100*i/length(sparseInds)))
    toc
end




mname=['/subMov_entropy_density_id' num2str(galID) '_' simDisplayName];
%mname2=[DEFAULT_MATFILE_DIR '/subMov_stars_id' num2str(galID) '_' simDisplayName];
save([DEFAULT_MATFILE_DIR '/' mname ],'F','F2','-v7.3')
%T_MATFILE_DIR '/' mname ],'F','F2','-v7.3')
fprintf('done')

galaxyMovieAssembly(F,mname);
galaxyMovieAssembly(F2,mname);

