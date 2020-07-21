%load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat');

cl=6;
new_env(cl);

prj='YZ';
ind=10;
bmap='*RdYlBu';
cc= brewermap(8,'Set1'); 

%shk=shockedge_a1;
%shkM1=shockedgeMask1_a1;
%shkM2=shockedgeMask2_a1;
%c1=shk{ind,4};
%c2=shk{ind,5};
%circ=cat(2,c1(shkM1(ind,:)),c2(shkM2(ind,:)));

%% boxx 8
boxx=8;
thick=0.1;
dil=7;

mkmap(boxx,'type','rho',prj,'clims',[8 14],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no',...
    'print','printag','papI');

mkmap(boxx,'type','f',prj,'clims',[-3 3],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no',...
    'print','printag','papI');

mkmap(boxx,'type','t',prj,'clims',[5 8],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no',...
    'print','printag','papI');

mkmap(boxx,'type','s',prj,'clims',[2 5],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no',...
    'print','printag','papI');

mkmap(boxx,'type','p',prj,'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no',...
    'print','printag','papI');


%% boxx 4 
boxx=4;
thick=0.05;
zoomb=[-1.5 -1.5 3];
tic=-1.5:0.5:1.5;
mkmap(boxx,'type','f',prj,'clims',[-3 3],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','zoom',zoomb,'marks',0.25.*get_rvir,...
    'xticks',tic,'yticks',tic,'print','printag','papI');

mkmap(boxx,'type','mach',prj,'clims',[0 2],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','zoom',zoomb,'marks',0.25.*get_rvir,...
    'xticks',tic,'yticks',tic,'print','printag','papI');


%% boxx 1 
boxx=1;
thick=0.025;

mkmap(boxx,'type','f',prj,'clims',[-3 3],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','marks',0.25.*get_rvir,...
    'print','printag','papI');

mkmap(boxx,'type','mach',prj,'clims',[0 2],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','marks',0.25.*get_rvir,...
    'print','printag','papI');

zoomb=[-0.35 -0.35 0.7];
arrow(1).start=[-0.26 0.3];
arrow(1).stop=[-0.115 0.25];
arrow(1).faceColor=cc(1,:);

arrow(2).start=[-0.27 0.115];
arrow(2).stop=[-0.15 0.18];
arrow(2).faceColor=cc(2,:); 
tic=-0.3:0.1:0.3;

mkmap(boxx,'type','t',prj,'clims',[7 8],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','zoom',zoomb,'marks',0.25.*get_rvir,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papI');

mkmap(boxx,'type','p',prj,'clims',[-1 1.5],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','zoom',zoomb,'marks',0.25.*get_rvir,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papI');


%% CL 14 
cl=14;
prj='YZ';
new_env(cl);

boxx=2;
mkmap(boxx,'type','f',prj,'clims',[-3 3],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','marks',0.25.*get_rvir,...
    'print','printag','papI');

boxx=1;
mkmap(boxx,'type','f',prj,'clims',[-3 3],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','marks',0.25.*get_rvir,...
    'print','printag','papI');


%% cluster 3 
cl=3;
prj='XY';
new_env(cl);
thick=0.05;
boxx=4;
zoomb=[-1.6 -1.6 3.2];
tic=-1.5:0.5:1.5;
mkmap(boxx,'type','f',prj,'clims',[-3 3],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','marks',0.25.*get_rvir,'zoom',zoomb,...
    'xticks',tic,'yticks',tic,'print','printag','papI');


new_env(cl,'a06');

boxx=4;
zoomb=[-0.85 -0.85 1.7];
tic=-1.5:0.5:1.5;
mkmap(boxx,'type','f',prj,'clims',[-3 3],'proper','vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','marks',0.25.*get_rvir,'zoom',zoomb,...
    'xticks',tic,'yticks',tic,'print','printag','papI');

%% cluster 11 
cl=11;
prj='XY';
new_env(cl);
thick=0.05;
boxx=2;

mkmap(boxx,'type','f',prj,'clims',[-3 3],'vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','marks',0.25.*get_rvir,...
    'print','printag','papI');


new_env(cl,'a06');

boxx=2;


mkmap(boxx,'type','f',prj,'clims',[-3 3],'proper','vfield','dilute',dil,'thick',thick,...
    'brewer',bmap,'labels','half','title','no','marks',0.25.*get_rvir,...
    'print','printag','papI');
