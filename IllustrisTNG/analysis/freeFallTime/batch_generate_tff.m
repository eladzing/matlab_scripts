%%  batch  - run the freefalltime calculation on all snapshots
% up to z=5. 

%% perliminaries - simulation environment must be set.


snapList=[17,21,25,33,40,50,59,67,72,78,84,91,99];

%% run over snapshots 

for i=1:length(snapList)
    
    snap=snapList(i);
    
    fprintf('Generating free-fall time profile for snapshot %i',snap);
    
    perlimFlag=true;
    
    freefalltime;
    
    clear tffProfile
    
end