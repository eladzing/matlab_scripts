%% recreate flux spheres. - obsolete - redid this while running batch_make_profile

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
units 
for i=list(1)
    for axp={'a1','a06'}
        type='csf';
        
        new_env(i,axp{1},type )
        global HALO_PATH
        global aexpn
        % calculate flux
        for box=[1 2 4 8];
           
           flux=new_flux(box); 
           vr=Vr_full(box);
           vr=cart2sphere(vr);  % in km/sec
           
           rog=RHOG_sphere(box); %in Msun/Mpc^3
           ds=ds_sphere(box); % in Mpc^2
           
           flux=rog.*vr.*ds.*(km/Mpc*yr); % in Msun/yr 
           
           save_cube(vr, HALO_PATH, sprintf('Vr_sphere_%s_%d',aexpn,box));
           save_cube(flux, HALO_PATH, sprintf('flux_sphere_%s_%d',aexpn,box));
           
            
        end
        
        
    end
end
