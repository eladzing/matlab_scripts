new_env('CL6','csf','a1')

global aexpn

path= '/home/eladzing/work/sshfs/sungate1/poteng/data/CL6/CSF';

for i=[1 2 4 8]
    
    [pdm pgs pdot ro]=poteng_temp('pot',i);
    epg=(pgs+pdm).*ro;
    edot=ro.*pdot;
    
    save_cube(cart2sphere(pdm), path, sprintf('potdm_sphere__%s_%d',aexpn,i));
    save_cube(cart2sphere(pgs), path, sprintf('potgas_sphere__%s_%d',aexpn,i));
    save_cube(cart2sphere(pdot), path, sprintf('potdot_sphere__%s_%d',aexpn,i));
    
    save_cube(cart2sphere(epg), path, sprintf('epgas_sphere_%s_%d',aexpn,i));
    save_cube(cart2sphere(edot), path, sprintf('epdot_sphere_%s_%d',aexpn,i));

    clear pgs pdot pdm ro epg edot                       
    
    
%     case 2
%         %if ~exist('hubflag','var')switch box
%         case 1
%             epg1=(pgs+pdm).*ro;
%             edot1=ro.*pdot;
%         case 2
%             epg2=(pgs+pdm).*ro;
%             edot2=ro.*pdot;
%         case 4
%             epg4=(pgs+pdm).*ro;
%             edot4=ro.*pdot;
%         case 8
%             epg8=(pgs+pdm).*ro;    %if ~exist('hubflag','var')    %if ~exi    %if ~exist('hubflag','var')st('hubflag','var')
%             edot8=ro.*pdot;
%     end

    
end