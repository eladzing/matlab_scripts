%% jf beyond r200 - spshaback or not
% how many of the JF found beyond r200 have not crossed r200 yet
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
%distMask15=rpos>=1.5;
%distMask2=rpos>=2;
vv=sqrt(Units.GG.*galProps.hostM200c./galProps.hostR200c);
vr=galProps.vrad'./vv';

mask=distMask & maskJF;

% %% plot position vs velocity
% scatter(rpos(distMask),vr(distMask),'.b')
% hold on
% scatter(rpos(distMask & maskJF),vr(distMask & maskJF),'.r')
% plot([1 10],[0 0],':k')

%%
treeType='SubLink_gal';
treeFields={'SnapNum','SubfindID','SubhaloPos',...
    'SubhaloID',...
    'GroupFirstSub','GroupPos','Group_R_Crit200'};
%'SubhaloSFRinRad','SubhaloMassInRadType','Group_M_Crit200',

rposOut=rpos(distMask);
rposMin=-1.*ones(size(rposOut));
timeLastInRv=rposMin;
timeLastInHalfRv=rposMin;
timeMin=rposMin;

indx=find(distMask);

%% find universal sound_crossing time 

zr=illustris.utils.snap2redshift(objectTable.snap(indx));
rhoc=rho_crit(zr,'cosmo',cosmoStruct);
tsc=sqrt((3/(4*pi*5/3)).*(Units.G.*200*rhoc'.*Units.Ms./Units.Mpc^3).^-1)./Units.Gyr;
clear rhoc;



%%
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
        
        %outskirt_histories.hist(i).tableIndx=ii;
        
       
                
        
        prc=floor(cnt./sum(mask)*100);
        if (prc>=stepNext)
            
            fprintf(' *** %i %% of histories completed  *** \n',prc);
            stepNext=stepNext+step;
        end
        
        tree = illustris.sublink.loadTree(bp,galProps.snap(ii),galProps.subfind(ii),treeFields,true,treeType);
        
        zr=illustris.utils.snap2redshift(tree.SnapNum);
        rposHist=findDistance(tree.SubhaloPos(:,:),tree.GroupPos(:,:),LBox,3)./tree.Group_R_Crit200'; %.*unitFactors.lengthUnit;
        
        isSat=tree.GroupFirstSub~=tree.SubfindID;
        
        zlast1=zr(find(rposHist(isSat)<1,1,'first'));
        zlast05=zr(find(rposHist(isSat)<0.5,1,'first'));
        
        if ~isempty(zlast1)
            timeLastInRv(i)=redshift2time(zr(1),'cosmo',cosmoStruct).age-...
                redshift2time(zlast1,'cosmo',cosmoStruct).age;  % time since last time it was within r200 of host.
        end
        
        if ~isempty(zlast05)
            timeLastInHalfRv(i)=redshift2time(zr(1),'cosmo',cosmoStruct).age-...
                redshift2time(zlast05,'cosmo',cosmoStruct).age;  % time since last time it was within r200 of host.
        end
        
%         outskirt_histories.hist(i).zr=zr;
%         outskirt_histories.hist(i).rpos=rposHist;
%         outskirt_histories.hist(i).isSat=isSat;
    end
end


%outskirt_histories.mask=mask;
outskirt_histories.tableIndices=indx;
outskirt_histories.isJF=maskJF(distMask);
outskirt_histories.tsc=tsc;
outskirt_histories.zred=zr;
outskirt_histories.rpos=rposOut;
outskirt_histories.timeLastInRv=timeLastInRv;
outskirt_histories.timeLastInHalfRv=timeLastInHalfRv;


%% save to file

fname=sprintf('outskirt_All_backsplash');

save([DEFAULT_MATFILE_DIR '/' fname],'outskirt_histories','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);

