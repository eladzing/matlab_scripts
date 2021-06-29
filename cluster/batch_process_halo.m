function batch_process_halo(halopath, haloformat, clustername)

global FILE_FORMAT;
FILE_FORMAT = sprintf('%s/%s', halopath, haloformat)
global FILE_FORMAT_SPHERE;
FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat')

%some cart2sphere preprocessing
for res = [8 4 2 1]
    try
        if (~exist(sprintf('%s/%s.mat', halopath, sprintf('RHOTOT_sphere_%d', res)), 'file'))
            disp(sprintf('cart2sphere, RHOTOT: %d', res))
            save_cube(cart2sphere(RHOTOT(res)), halopath, sprintf('RHOTOT_sphere_%d',res));
        end
        if (~exist(sprintf('%s/%s.mat', halopath, sprintf('RHOG_sphere_%d', res)), 'file'))
            disp(sprintf('cart2sphere, RHOG: %d', res))
            save_cube(cart2sphere(RHOG(res)), halopath, sprintf('RHOG_sphere_%d',res));
        end
        if (~exist(sprintf('%s/%s.mat', halopath, sprintf('T_sphere_%d', res)), 'file'))
            disp(sprintf('cart2sphere, T: %d', res))
            save_cube(cart2sphere(T(res)), halopath, sprintf('T_sphere_%d',res));
        end
        if (~exist(sprintf('%s/%s.mat', halopath, sprintf('Vr_sphere_%d', res)), 'file'))
            disp(sprintf('cart2sphere, Vr: %d', res))
            save_cube(cart2sphere(Vr(res)), halopath, sprintf('Vr_sphere_%d',res));
        end
    catch
        disp(lasterror)
    end
end

%profiles
disp(sprintf('profiles: %d', 8))
found_rvir = batch_profiles(8, halopath, clustername, 0);
if (found_rvir)
    for res = [4 2 1]
        profile_filename = sprintf('%s/virial%d.mat',halopath, res);
        if (~exist(profile_filename, 'file'))
            disp(sprintf('profiles: %d', res))
            try
                batch_profiles(res, halopath, clustername, 1);
            catch
                disp(lasterror)
            end
        end
    end
end
close all

%flux
for res = [8 4 2 1]
    try
        flux_filename = sprintf('%s/fluxes%d.mat',halopath, res);
        if (~exist(flux_filename, 'file'))
            disp(sprintf('flux: %d', res))
            batch_flux(res, halopath, clustername);
        end
    catch
        disp(lasterror)
    end
end
close all

%bird
for res = [8 4 2 1]
    try
        bird_filename = sprintf('%s/bird%d.mat',halopath, res);
        if (~exist(bird_filename, 'file'))
            disp(sprintf('bird: %d', res))
            [result_radii result_radii_mass result_radii_cnt result_map result_vrho result_vtemp result_count] = work_bird(res);
            save(bird_filename, 'result_radii', 'result_radii_mass', 'result_radii_cnt', 'result_map', 'result_vrho', 'result_vtemp', 'result_count');
        end
    catch
        disp(lasterror)
    end
end

%turbulence
try
    turb_filename = sprintf('%s/turbulence.mat',halopath);
    if (~exist(turb_filename, 'file'))
        disp('turbulence')
        batch_turbulence(halopath, clustername);
    end
catch
    disp(lasterror)
end

close all
