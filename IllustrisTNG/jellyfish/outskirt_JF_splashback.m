%% jf beyond r200 - spshaback or not
% how man y of the JF found beyond r200 have not crossed r200 yet
% and how many are 'splashback'


%% load object tables

global DEFAULT_MATFILE_DIR
% load object table
load([ DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat']);
load([DEFAULT_MATFILE_DIR '/jf_galProperties_CJF.mat']);

global cosmoStruct
global LBox
%%
units;
maskJF=objectTable.scoreWeighted>=0.8;


rpos= (galProps.rpos./galProps.hostR200c)';
distMask=rpos>=1;
distMask15=rpos>=1.5;
distMask2=rpos>=2;
vv=sqrt(Units.GG.*galProps.hostM200c./galProps.hostR200c);
vr=galProps.vrad'./vv';

mask=distMask & maskJF;

%% plot position vs velocity
scatter(rpos(distMask),vr(distMask),'.b')
hold on
scatter(rpos(distMask & maskJF),vr(distMask & maskJF),'.r')
plot([1 10],[0 0],':k')

%%
treeType='SubLink_gal';
treeFields={'SnapNum','SubfindID','SubhaloPos',...
    'SubhaloID',...
    'GroupFirstSub','GroupPos','Group_R_Crit200'};
%'SubhaloSFRinRad','SubhaloMassInRadType','Group_M_Crit200',

rposJF=rpos(mask);
rposMin=-1.*ones(size(rposJF));
timeLastInRv=rposMin;
timeMin=rposMin;

indx=find(mask);

sims=unique(galProps.sim);

cnt=0;
step=5;
stepNext=0;
for k=1:2
    
    bp=illustris.set_env(sims(k));
    
    for i=1:length(indx)
        ii=indx(i);
        if galProps.sim(ii)~=sims(k)
            continue;
        end
        
        cnt=cnt+1;
        
        res.hist(i).tableIndx=ii;
        
        
        prc=floor(cnt./sum(mask)*100);
        if (prc>=stepNext)
            
            fprintf(' *** %i %% of histories completed  *** \n',prc);
            stepNext=stepNext+step;
        end
        
        tree = illustris.sublink.loadTree(bp,galProps.snap(ii),galProps.subfind(ii),treeFields,true,treeType);
        
        zr=illustris.utils.snap2redshift(tree.SnapNum);
        rposHist=findDistance(tree.SubhaloPos(:,:),tree.GroupPos(:,:),LBox,3)./tree.Group_R_Crit200'; %.*unitFactors.lengthUnit;
        
        isSat=tree.GroupFirstSub~=tree.SubfindID;
        
        [rposMin(i),ix]=min(rposHist(isSat)); % minimal distance form host when sat
        
        timeMin(i)=redshift2time(zr(1),'cosmo',cosmoStruct).age-...
            redshift2time(zr(ix),'cosmo',cosmoStruct).age;    % time of minimal distance 
        
        zlast=zr(find(rposHist<1,1,'first'));
        if ~isempty(zlast)
        timeLastInRv(ii)=redshift2time(zr(1),'cosmo',cosmoStruct).age-...
            redshift2time(zlast,'cosmo',cosmoStruct).age;  % time since last time it was within r200 of host.
        end
        
        res.hist(i).zr=zr;
        res.hist(i).rpos=rposHist;
        res.hist(i).isSat=isSat;
    end
end


res.mask=mask;
res.indices=indx;
res.rposJF=rposJF;
res.rposMin=rposMin;
res.timeMin=timeMin;
res.timeLastInRv=timeLastInRv;


%% save to file

fname=sprintf('outskirt_JF_backsplash');

save([DEFAULT_MATFILE_DIR '/' fname],'centralHist','done','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);

