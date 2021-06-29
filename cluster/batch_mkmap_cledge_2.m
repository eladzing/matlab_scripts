

list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
mask=[ 0   0   1   0   0   1   0   0  0  1  0  0  1  1  1  0];
%mask=[ 0   0   0   0   0   0   0   0  0  1  0  0  0  1  1  0];
%mask=true(size(list));
%mask(10)=true;

%list=[103 106 6 10 11 14];
%mask=[ 1   1  1  1  1  1];
%boxx=8;

%projs=ones(length(list),3);

% 1 - yz 2- zx 3- xy
 projs(3,:)=[0 1 0]; %103
 projs(6,:)=[0 0 1]; %106
 projs(10,:)=[1 1 1]; %6
 projs(13,:)=[0 1 0]; %10
 projs(14,:)=[0 1 0]; %11
 projs(15,:)=[0 0 1]; %14


boxx=8;
% z=0
for i=1:length(list);
    %continue % THIS IS HERE SINCE THESE PICS ARE DONE
    if~mask(i)
        continue;
    end
    for j=1:2
        switch j
            case 1
                new_env(sprintf('CL%d',list(i)),'csf','a1');
            case 2
                new_env(sprintf('CL%d',list(i)),'csf','a06');
        end
        
        global DEFUALT_PRINTOUT_DIR;

        printoutdir=sprintf('%s/%s',DEFUALT_PRINTOUT_DIR,'rps_full');

        rog=RHOG(boxx);
        rpscube=rog.*get_vvir().^2;
        vx = Vx(boxx);
        vy = Vy(boxx);
        vz = Vz(boxx);
        
        mkmap_rps(boxx,'v','verbose','rpscube',rpscube,'rog',rog,...
            'vxcube',vx,'vycube',vy,'vzcube',vz,'proj',projs(i,:),...
            ...%'contour','mvir','cont_val',1e12,...
            'type','mv_strip_inv','strip_model','mv_strip','strip_param',1e11,...
            'vfield','dilute',8,'clims',[0 1],...
            'printag','delete','proper','printout',printoutdir)
        
%          mkmap_rps('rpscube',rpscube,'rog',rog,...
%             'vxcube',vx,'vycube',vy,'vzcube',vz,'proj',projs(i,:),...
%             'title','RPS of Exp. Disk','boxx',boxx,'log',...%'contour','mvir','cont_val',1e12,...
%             'type','expdisk_sig_5090','strip_model','expdisk_sigma','strip_param',[4e7 0.1 1.0],...
%             'rps_v','vvir','vfield','yes','dilute',8,...
%             'printag','rps_rs12','proper','printout',printoutdir)
%         
%         
       % mkmap_rps('data',rpscube,'rpscube',rpscube,'rog',rog,'vxcube',vx,'vycube',vy,'vzcube',vz,'proj',projs(i,:),'title','RPS','boxx',boxx,'log','on','contour','rstrip','cont_val',0.1,'rps_v','vvir','vfield','yes','dilute',8,'print','yes','printag','rps_mv10','proper','yes','printout',printoutdir)
    end
   
    %close all
    
end

% 
% % z=0.6
% for i=1:length(list);
%     %continue % THIS IS HERE SINCE THESE PICS ARE DONE
%     if~mask(i)
%         continue;
%     end
%     
%     
%     
%     mkmap_rps('type','rps','proj',projs(i,:),'title','RPS','boxx',boxx,'log','on','contour','mvir','cont_val',1e12,'rps_v','vvir','vfield','yes','dilute',8,'print','yes','printag','rps_mv','proper','yes')
%     mkmap_rps('type','rps','proj',projs(i,:),'title','RPS','boxx',boxx,'log','on','contour','rstrip','cont_val',0.1,'rps_v','vvir','vfield','yes','dilute',8,'print','yes','printag','rps_mv','proper','yes')
%     close all
% end

%
%
% for i=1:length(list);
%     if~mask(i)
%         continue;
%     end
%
%    new_env(sprintf('CL%d',list(i)),'csf','a06');
%
%    for boxx=[8];
%           if (i==10)
%             lf='on';
%             else
%             lf='half';
%           end
%           circs=edge_a06(i,:).*rv_a06(i);
%           mkmap('type','Entropy','box',boxx,'clims',[-3 3],'proj',projs(i,:),'thick',1/256*6,'marks',circs,'title','','labels',lf,'print','yes','printag','ent_a06_cledge','proper','yes')
%           mkmap('type','temperature','box',boxx,'clims',[4 8],'proj',projs(i,:),'thick',1/256*6,'marks',circs,'title','','labels',lf,'print','yes','printag','tmp_a06_cledge','proper','yes')
%           %mkmap('Entropy',boxx,[-3 3],projs(i,:),1/256*6,'novel',-7,circs,[],'','print','ent_a06_cledge')
%           %mkmap('Temperature',boxx,[4 8],projs(i,:),1/256*6,'novel',-7,circs,[],'','print','tmp_a06_cledge')
% %          circs=edge2(i,:);
% %          mkmap('Entropy',boxx,[-3 3],[1 1 1],1/256*6,'novel',-7,circs,[],'Edge by $max(S)$','print','ent_a1_testedge2',printoutdir)
% %          mkmap('Temperature',boxx,[4 8],[1 1 1],1/256*6,'novel',-7,circs,[],'Edge by $max(S)$','print','tmp_a1_testedge2',printoutdir)
%
%        %mkmap('rps',8,[15 20],[1 1 1],boxx/256*6,'no',-4,[],1,'$z=0$','print','rps_a1')
%    end
%
%    %close all
% end
