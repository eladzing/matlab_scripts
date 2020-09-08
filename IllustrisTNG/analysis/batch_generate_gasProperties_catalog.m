%%  batch  - run the generate_gasProperties 
% up to z=5. 

%% perliminaries - simulation name must be set (simName)

%bp=illustris.set_env(simName);
snapList=[17,21,25,33,40,50,59,67,72,78,84,91,99];

%% run over snapshots 
if ~exist('startInd','var')
    startInd=1;
end
for i=startInd:length(snapList)
    
    snap=snapList(i);
    
    fprintf('Generating gasProperties catalog for snapshot %i \n',snap);
    
    readFlag=true;
    
    generate_gasProperties_forCatalog
    
    clear PropStruct 
    
end