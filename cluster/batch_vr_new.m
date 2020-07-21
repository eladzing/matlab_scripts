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
        vr_sph=cart2sphere(Vr_full(boxx));
        ro_sph=RHOG_sphere(boxx);
        
        flux=new_flux(boxx,'vr_sph',vr_sph,'ro_sph',ro_sph);
        
        %save_cube(vr_sph, HALO_PATH, sprintf('Vr_sphere_%s_%d',aexpn,boxx));
        save_cube(flux, HALO_PATH, sprintf('flux_sphere_%s_%d',aexpn,boxx));
    end
end