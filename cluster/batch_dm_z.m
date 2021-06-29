%function [md md_in md_out md_in_strong vrp vrp_in vrp_out vrp_in_strong r_prof]= flux_profs_new(varargin)
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
mask=true(length(list));

for i=1:length(list)
    c1=clock;
    if~mask(i)
        continue
    end
    for k=2
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
            mkmap('box',boxx,'marks',get_rvir,'type','flux','proj',[1 1 1],...
                'clims',[-3 3],'vfield','yes','proper','yes',...
                'print','yes','format','both','printag','flux',...
                'printout',sprintf('%s/flux_maps',DEFAULT_PRINTOUT_DIR))
            
            
            
         close all   
        end
    end
    c2=clock-c1
end