%% draw flux from outside in, in sphericl shells

global NCELL

ni=NCELL*5/2;


cu=flux_sphere(boxx);
[mesh_r mesh_phi mesh_theta] = sphere_grid(boxx);
read=false;
i=ni-1;
step=60;
bulk=NCELL/2*4;
while i>1boxx=8;

    
    rad=mesh_r(i-bulk+128,1,1)*0.7
    
    % move to smaller box
    if rad<boxx/4
        boxx=0.5*boxx;
        % read smaller box
        cu=flux_sphere(boxx);
        [mesh_r mesh_phi mesh_theta] = sphere_grid(boxx);
        
        switch boxx
            case 4
                bulk=NCELL/2*3;
            case 2
                bulk=NCELL/2*2;
            case 1
                bulk=0;
        end
        
    end
    
    ind=i-bulk+128;
    
    
    shel=squeeze(cu(ind,:,:));
    sheli=shel;
    sheli(shel>0)=1e-30;
    mphi=squeeze(mesh_phi(ind,:,:));
    mthet=squeeze(mesh_theta(ind,:,:));
    
    plot_hammer_shell(log10(abs(sheli)),mthet,mphi,400)
    caxis([-2 1])
     pause
    i=i-step;
end