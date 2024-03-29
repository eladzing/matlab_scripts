%function [md md_in md_out md_in_strong vrp vrp_in vrp_out vrp_in_strong r_prof]= flux_profs_new(varargin)
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
mask=true(length(list));

for i=1:length(list)
    if~mask(i)
        continue
    end
    new_env(sprintf('CL%d',list(i)),'csf','a1','win');
    
    
    
    global HALO_PATH
    global aexpn
   
    
    
    for j=1:4;
        boxx=2^(j-1);
        %load(sprintf('%s/virial%d_%s', HALO_PATH, boxx,aexpn));
        %% find vr and flux
        
        
        
        s_sph=cart2sphere(S(boxx));
                
        %save_cube(vr_sph, HALO_PATH, sprintf('Vr_sphere_%s_%d',aexpn,boxx));
        save_cube(s_sph, HALO_PATH, sprintf('S_sphere_%s_%d',aexpn,boxx));
    end
end