snaps=[33    40    50    67    84    99];

global DEFAULT_MATFILE_DIR
global simDisplayName
global illUnits

for i= 5 %1:length(snaps)
    
    snap=snaps(i);
    fprintf('snap %i \n',snap);
    fprintf('reading fofs \n');
    
    fofs=illustris.groupcat.loadHalos(bp,snap);
    fprintf('reading subs \n');
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
    illustris.utils.set_illUnits(snap)
    
    
    stellar_mass=double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit);
    host_mass=double(fofs.Group_M_Crit200.*illUnits.massUnit);
    
    switch snap
        case 33
            
            for j=1:length(subfindList_33)
                mass33(j,1)=stellar_mass(subfindList_33(j,1)+1);
                mass33(j,2)=host_mass(subfindList_33(j,2)+1);
            end
        
            
            case 40
            
            for j=1:length(subfindList_40)
                mass40(j,1)=stellar_mass(subfindList_40(j,1)+1);
                mass40(j,2)=host_mass(subfindList_40(j,2)+1);
            end
        
            
            case 50
            
            for j=1:length(subfindList_50)
                mass50(j,1)=stellar_mass(subfindList_50(j,1)+1);
                mass50(j,2)=host_mass(subfindList_50(j,2)+1);
            end
        
            
            case 67
            
            for j=1:length(subfindList_67)
                mass67(j,1)=stellar_mass(subfindList_67(j,1)+1);
                mass67(j,2)=host_mass(subfindList_67(j,2)+1);
            end
        
            
            case 84
            
            for j=1:length(subfindList_84)
                mass84(j,1)=stellar_mass(subfindList_84(j,1)+1);
                mass84(j,2)=host_mass(subfindList_84(j,2)+1);
            end
        
            
            case 99
            
            for j=1:length(subfindList_99)
                mass99(j,1)=stellar_mass(subfindList_99(j,1)+1);
                mass99(j,2)=host_mass(subfindList_99(j,2)+1);
            end
    end
    
end

    fname=['jflist_masses_' simDisplayName '.mat'];
    fprintf('writing to file: %s  \n',fname);
    save([DEFAULT_MATFILE_DIR '/' fname],'mass33','mass40',...
        'mass50','mass67',...
        'mass84','mass99','-v7.3')
