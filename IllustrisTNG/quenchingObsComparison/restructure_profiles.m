%% short script to load profiles and assign them to data structures
global DEFAULT_MATFILE_DIR
list={'100','100-2','100-3',...
    '300','300-2','300-3',...
    '50','50-2','50-3','50-4'};


for i=1:length(list)
    
    fname=[DEFAULT_MATFILE_DIR '/ssfr_stellarMass_radProfiles_TNG' list{i} '.mat'];
    load(fname);
    
    switch list{i}
        case '100'
            TNG100(1)=profStruct;
            %             TNG100(1).ssfrAvg=ssfrAvg;
            %             TNG100(1).ssfrAvgC=ssfrAvgC;
            %             TNG100(1).starMass=starMass;
            %             TNG100(1).starMassC=starMassC;
            
        case '100-2'
            TNG100(2)=profStruct;
            
        case '100-3'
            TNG100(3)=profStruct;
            
        case '300'
            TNG300(1)=profStruct;
            
        case '300-2'
            TNG300(2)=profStruct;
            
        case '300-3'
            TNG300(3)=profStruct;
            
        case '50'
            TNG50(1)=profStruct;
            
        case '50-2'
            TNG50(2)=profStruct;
            
            
        case '50-3'
            TNG50(3)=profStruct;
            
            
        case '50-4'
            TNG50(4)=profStruct;
            
            
    end
    
    clear profStruct
end