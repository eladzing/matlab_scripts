cl=103;

%aa='a1';
%prj=[0 1 0]; %CL103
%prj=[1 0 0]; %CL6

%load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat')
load('/home/zinger/workProjects/matlab_scripts/cluster/mat_files/shockedge.mat')

switch cl
    case 103
        ind=3;prj=[0 1 0];
    case 6
        ind=10;prj=[1 0 0];
end


boxx=8;

for j=1:2
    
    switch j
        case 1
            aa='a1';
            shk=shockedge_a1;
            shkM1=shockedgeMask1_a1;
            shkM2=shockedgeMask2_a1;
        case 2
            aa='a06';
            shk=shockedge_a06;
            shkM1=shockedgeMask1_a06;
            shkM2=shockedgeMask2_a06;
    end
    new_env(cl,aa)
    rog=RHOG(boxx);
    rpscube=rog.*get_vvir().^2;
    vx = Vx(boxx);
    vy = Vy(boxx);
    vz = Vz(boxx);
    
    
    c1=shk{ind,4};
    c2=shk{ind,5};
    
    circ=cat(2,c1(shkM1(ind,:)),c2(shkM2(ind,:)));
    circ2=cat(2,c1,c2);
    
    
%     mkmap_rps_nfw(boxx,'proj',prj,'type','mv_strip_inv',...
%         'strip_model','mv_strip','strip_param',1e11,...
%         'rpscube',rpscube,'rog',rog,...
%         'vxcube',vx,'vycube',vy,'vzcube',vz,'marks',circ,...
%         'brewer','*RdYlBu','proper','vfield','dilute',8,...
%         'alfa',0.5,'clims',[0 100],'labels','half','title','no',...
%         'print','printag','nfw_papII2');
    
     mkmap(boxx,'xz','type','rps',...
        'rog',rog,'velcolor','k','clims',[-3 1.5],...
        'vxcube',vx,'vycube',vy,'vzcube',vz,'marks',circ,...
        'brewer','YlOrRd','proper','vfield','dilute',10,...
        'alfa',0.5,'labels','half','title','no',...
        'print','printag','rps_papII2','streamlines','off');
    
     
     
    % 8,'type','rps','alpha',0.5,'xz')
    
    
%     mkmap_rps_nfw(boxx,'proj',prj,'type','rstrip',...
%         'strip_model','rstrip','strip_param',0.1,...
%         'rpscube',rpscube,'rog',rog,...
%         'vxcube',vx,'vycube',vy,'vzcube',vz,'marks',circ,...
%         'brewer','*RdYlBu','proper','vfield','dilute',8,...
%         'clims',[8 12],'alfa',0.5,'labels','half','title','no',...
%         'print','printag','papIIDel');
end