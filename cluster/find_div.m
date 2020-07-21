%% calculate divergence of velocity field

typ='csf';aexp='a06';

list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 24];

for id=1:length(list)

   clname=sprintf('CL%d',list(id));
   new_env(clname,typ,aexp)
   
   
   %%find Vcm
    global hub;
    global HALO_PATH;
    global NCELL;
    
    load(sprintf('%s/virial%d_%s', HALO_PATH, 8,aexp));
    h=hub;
    rvfac=1.0;
    rvcm=RVIR*rvfac;
    vbox=ceil(2^ceil(log2(2*rvcm*h)));
    
    rc=mk_rcube(vbox,ones(NCELL,NCELL,NCELL));
    rv_ind=find(rc<=rvcm);
    [VcmX VcmY VcmZ] = V_Vcm_r(vbox,rv_ind,'hub'); 
    clear rv_ind rc vbox;
    for i=1:4
        boxx=2.^(i-1)
        
        [hubX hubY hubZ] = hubble_flow(boxx);
        Vxx = Vx(boxx)+hubX;
        Vyy = Vy(boxx)+hubY;
        Vzz = Vz(boxx)+hubZ;
    
        Vxx = Vxx - VcmX;
        Vyy = Vyy - VcmY;
        Vzz = Vzz - VcmZ;
    
        [meshY, meshX, meshZ] = meshgrid(1:NCELL, 1:NCELL, 1:NCELL);

        %convert to center origin coordinates
        meshX = meshX - (NCELL+1)/2;% -cm(1);
        meshY = meshY - (NCELL+1)/2;% -cm(2);
        meshZ = meshZ - (NCELL+1)/2;% -cm(3);
        % Fix Units (to be in Mpc)
        meshX = meshX * ((boxx/h)/NCELL);
        meshY = meshY * ((boxx/h)/NCELL);
        meshZ = meshZ * ((boxx/h)/NCELL);
    
        div=divergence(meshY,meshX,meshZ,Vxx,Vyy,Vzz);
    
        save_cube(div, HALO_PATH, sprintf('divV_%s_%d',aexp,boxx));
        clear div meshX meshY meshZ Vxx Vyy Vzz
    end
end
