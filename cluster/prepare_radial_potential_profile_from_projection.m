%% generate radial profiles of the potential projection 

new_env(107);

nbins=200;
%% read projection 
for i=1:3
    switch i
        case 1
            projName='cl107_pot_n500_projXY.dat';
        case 2
            projName='cl107_pot_n500_projXZ.dat';
        case 3
            projName='cl107_pot_n500_projYZ.dat';
    end
    res = read_pot_proj( projName);
    
%% generate profile
prof(i) = mk_radial_profile_2d(res.arr,'bottomLeft',res.bottomLeft,'side',side,'nbins',nbins);

%% write to file

end