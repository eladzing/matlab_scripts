function outMask = get_fof_icm(fofs,id,snap)
%GET_FOF_ICM find the indices of the fof which correspond to the main subhalo
%and "fuzz"
%   For a given fof, get a list of gas cells which are either part of the
%   main subhalo or part of the FOF "fuzz", i.e., are not connected with
%   any smaller subhalo. In essences, this means excising all the subhalos
%   which are not the main subhalo

%% get cell indices of all FOF cells

global BASEPATH
groupOffset = illustris.snapshot.getSnapOffsets(BASEPATH,snap,id,'Group');% get the relevant particle indices

if groupOffset.lenType(1)
    fprintf('%s - FoF ID %i has no gas elements \n',current_function().upper,id);
    return
end

groupFirstInd=groupOffset.offsetType(1)+1;
%groupLastInd=groupFirstInd+int64(groupOffset.lenType(1))-1;

outMask=true(1,groupOffset.lenType(1));

%% test
fofGas=illustris.snapshot.loadHalo(BASEPATH, snap, id, 'gas',{'ParticleIDs'});

%% generate a list of subhalos
subIDs=fofs.GroupFirstSub(id+1) + (0:fofs.GroupNsubs(id+1)-1);

%% loop over list except for main sub
for i=2:length(subIDs)
    
    
    subOffset= illustris.snapshot.getSnapOffsets(BASEPATH,snap,subIDs(i),'Subhalo');% get the relevant particle indices
    if subOffset.lenType(1)==0 % skip if subhalo has no gas 
        continue
    end
    
    subFirstInd=subOffset.offsetType(1)+1;
    subLastInd=subFirstInd+int64(subOffset.lenType(1))-1;
    
    inds=subFirstInd:subLastInd ;
    inds=inds - groupFirstInd +1 ;
    
    outMask(inds)=false;
    
    %% test
    try
        satGas=illustris.snapshot.loadSubhalo(BASEPATH, snap, subIDs(i), 'gas',{'ParticleIDs'});
        if any(satGas~=fofGas(inds))
            error('particle ids do not match up')
        end
    catch
        error('subhalo ID %i,  %i \n',i,subIDs(i))
    end
end


end

