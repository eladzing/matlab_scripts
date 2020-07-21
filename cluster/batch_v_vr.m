%function [md md_in md_out md_in_strong vrp vrp_in vrp_out vrp_in_strong r_prof]= flux_profs_new(varargin)
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
mask=false(length(list));
mask(10)=true;

for i=1:length(list)
    c1=clock;
    if~mask(i)
        continue
    end
    for k=1
        if k==1
            new_env(sprintf('CL%d',list(i)),'csf','a1','win');
        else
            new_env(sprintf('CL%d',list(i)),'csf','a06','win');
        end
        
        %global HALO_PATH
        %global aexpn
        global DEFAULT_PRINTOUT_DIR
        
        
        
        for j=1:4;
            boxx=2^(j-1);
            
            if read
                rog=RHOG(boxx);
                vx = Vx(boxx);
                vy = Vy(boxx);
                vz = Vz(boxx);
            end
            
            mkmap('box',boxx,'marks',get_rvir,'type','velocity','proj',[1 1 1],...
                'vxcube',vx,'vycube',vy,'vzcube',vz,'norm',get_vvir,...
                'vfield','yes','proper','yes','log','yes',...
                'print','yes','format','both','printag','vtot_norm',...
                'printout',sprintf('%s/vel_maps',DEFAULT_PRINTOUT_DIR))
            
            mkmap('box',boxx,'marks',get_rvir,'type','vr','proj',[1 1 1],...
                'vxcube',vx,'vycube',vy,'vzcube',vz,'norm',get_vvir,...
                'vfield','yes','proper','yes','log','no',...
                'print','yes','format','both','printag','vr_norm',...
                'printout',sprintf('%s/vel_maps',DEFAULT_PRINTOUT_DIR))
            
            
            
            close all
        end
    end
    c2=clock-c1
end