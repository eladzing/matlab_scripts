function [outputArg1,outputArg2] = get_fof_icm(id,snap)
%GET_FOF_ICM find the indices of the fof which correspond to the main subhalo
%and "fuzz"
%   For a given fof, get a list of gas cells which are either part of the
%   main subhalo or part of the FOF "fuzz", i.e., are not connected with
%   any smaller subhalo. In essences, this means excising all the subhalos
%   which are not the main subhalo 

%% get cell indices of all FOF cells

global BASEPATH
fofCells = illustris.snapshot.getSnapOffsets(BASEPATH,snap,id,'Group');% get the relevant particle indices

%% generate a list of subhalos 
subIDs=fofs.GroupFirstSub(id+1) + (0:fofs.GroupNsubs(id+1)-1);

%% loop over list 
for i=1:length(subIDs)
           

    subCells= illustris.snapshot.getSnapOffsets(BASEPATH,snap,subIDs(i),'Subhalo');% get the relevant particle indices

% if subhalo is main subhalo include the particles 

% if part of satellite subhalo, remove it 



outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

