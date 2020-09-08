%%  batch  - run the generate_gasProperties 
% up to z=5. 

%% perliminaries - simulation name must be set (simName)

%bp=illustris.set_env(simName);
snapList=[17,21,25,33,40,50,59,67,72,78,84,91,99];

%% run over snapshots 
if ~exist('startInd','var')
    startInd=1;
end

global DEFAULT_MATFILE_DIR
global simDisplayName
compNames=["Gal", "CGMin", "CGMout", "CGMall" "Sub"];

for i=startInd:length(snapList)
    
    snap=snapList(i);
    
    fprintf('Reading in gasProperties mat file for snapshot %i \n',snap);
    
    fname=sprintf('gasProperties_snp%s_%s.mat',num2str(snap),simDisplayName);
    load([DEFAULT_MATFILE_DIR '/' fname],'PropStruct')

    for fld=compNames
        
        outStruct=PropStruct.(fld);
                
        catName=sprintf('Subhalos_%s_PhysicalGasProperties',fld);
        
        folder=['PhysicalGasProperties/' simDisplayName];
        
        illustris.utils.write_catalog(outStruct,snap,'name',catName,...
            'path','default','folder',folder,'v');
    
    end
    clear PropStruct 
    
end