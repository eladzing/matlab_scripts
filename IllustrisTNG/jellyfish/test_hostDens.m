%% script to examine

%% load data
global illUnits

global DEFAULT_MATFILE_DIR
% load object table
load([ DEFAULT_MATFILE_DIR '/cosmic_jellyfish_objectTable.mat']);

% load galaxy properties
load([DEFAULT_MATFILE_DIR '/jf_galProperties_CJF.mat']);

%% perliminaries
sims={'TNG50','TNG100'};
%maskJF=objectTable.scoreWeighted>=0.8;

snaps=unique(objectTable.snap);

mv=galProps.hostM200c;
rv=galProps.hostR200c;
rp=galProps.rpos./rv;

objList=find(rp>0);
simList=galProps.sim(objList)';
snapList=galProps.snap(objList)';
zred=illustris.utils.snap2redshift(snapList);
% prepare output
res.objectList=objList;
res.snap=snapList;
res.sim=simList;
res.hostID=galProps.hostID(objList);
res.hostM200c=mv(objList);
res.hostR200c=rv(objList);
res.hostRpos=galProps.rpos(objList);

cv0=cvir_Mvir_200(mv(objList),zred,'200');

res.hostDens=nfw.nfwRho(rp(objList),cv0); % not normalized! 


fprintf('Examining %i objects \n',length(objList));

for k=1:length(sims)
    bp=illustris.set_env(sims{k});
    global LBox
    
    for i=length(snaps)
        snp=snaps(i);
        zr=illustris.utils.snap2redshift(snp);
        listInd=find(simList==sims{k} & snapList==snp); % indices of relevant objects for the sim/snap in objList
        subList=objList(listInd); 
        
         fprintf('%s , %i) snap %i ',sims{k},i,snp);
        
        % if no objects are relavent go to next iteration
        if isempty(subList)
            fprintf(' .... no objects \n');
            continue
        else
            fprintf(' ... %i objects \n',length(subList));
        end
        
        % load gal/host catalogs
        [subs,fofs,subsInfo]=illustris.loadFofSub(snp);
        illustris.utils.set_illUnits(snp);
        
        
        mm=double(fofs.Group_M_Crit200).*illUnits.massUnit;
        rr=double(fofs.Group_R_Crit200).*illUnits.lengthUnit;
        gasFrac=double(fofs.GroupMassType(illustris.partTypeNum('gas')+1,:))./...
            double(fofs.GroupMassType(illustris.partTypeNum('dm')+1,:));
        
        msk=mm>=1e10;
                
        iii=find(msk);
        mmm=mm(msk);
        rrr=rr(msk);
        gf=gasFrac(msk);
        hostPos=double(fofs.GroupPos(:,msk));
        cv=cvir_Mvir_200(mmm,zr,'200');
        %rhoHalo=mmm./(4*pi/3*rrr.^3).*illUnits.densityUnit;
        
        
        for j=1:length(subList)
            ind=subList(j);
            sid=objectTable.subfind(ind);
            hid=galProps.hostID(ind);
           
            galpos=double(subs.SubhaloPos(:,sid+1));
            
            
            rpos=findDistance(galpos,hostPos,LBox,3).*illUnits.lengthUnit./rrr;
            
            dens=nfw.nfwRho(rpos,cv);
            densG=nfw.nfwRho(rpos,cv).*gf;
            
            
            %% based on NFW alone 
            [dmax,ix]=max(dens);
            maxInd=iii(ix); % index of maximal host in FoF catalog
            res.closeID(listInd(j))=maxInd-1;
            res.closeM200c(listInd(j))=illUnits.massUnit.*...
                double(fofs.Group_M_Crit200(maxInd));
            res.closeR200c(listInd(j))=illUnits.lengthUnit.*...
                double(fofs.Group_R_Crit200(maxInd));
            res.closeRpos(listInd(j))=rpos(ix);
            res.closeDens(listInd(j))=dmax;
            res.closeHostFlag(listInd(j))=(maxInd-1)==hid;
            
               %% based on NFW and global mass fraction
            [dmaxG,ixG]=max(densG);
            maxIndG=iii(ixG); % index of maximal host in FoF catalog
            res.closeID_gf(listInd(j))=maxIndG-1;
            res.closeM200c_gf(listInd(j))=illUnits.massUnit.*...
                double(fofs.Group_M_Crit200(maxIndG));
            res.closeR200c_gf(listInd(j))=illUnits.lengthUnit.*...
                double(fofs.Group_R_Crit200(maxIndG));
            res.closeRpos_gf(listInd(j))=rpos(ixG);
            res.closeDens_gf(listInd(j))=dmaxG;
            res.closeHostFlag_gf(listInd(j))=(maxIndG-1)==hid;
            
            
            

            
        end
        
    end
end

compareDens=res;

%% write to file 
fname=sprintf('jf_compareDEnsity_NFW_gasFrac_CJF.mat');
save([DEFAULT_MATFILE_DIR '/' fname],'compareDens','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);


