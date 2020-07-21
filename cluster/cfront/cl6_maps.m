load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat');

cl=6;
new_env(cl);

prj='YZ';
ind=10;
bmap='*RdYlBu';


shk=shockedge_a1;
shkM1=shockedgeMask1_a1;
shkM2=shockedgeMask2_a1;
c1=shk{ind,4};
c2=shk{ind,5};
circ=cat(2,c1(shkM1(ind,:)),c2(shkM2(ind,:)));

% boxx 8
thick=0.5;

mkmap(8,'type','rho',prj,'clims',[7 14],'vfield','dilute',6,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

mkmap(8,'type','numberdensity',prj,'clims',[-9 -2],'vfield','dilute',6,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

mkmap(8,'type','dm',prj,'clims',[7 15],'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

%%boxx 4 
thick=0.100;
mkmap(4,'type','f',prj,'clims',[-3 3],'vfield','dilute',6,'thick',thick,...
    'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

%%boxx 2
bmap='*Spectral';
thick=0.05;
mkmap(2,'type','vel',prj,'clims',[1.8 3.2],'log','vfield','dilute',6,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'bartag','$\log(v)\,[\mathrm{km/sec}]$','print','printag','papIII');

mkmap(2,'type','Mach',prj,'clims',[0 2],'vfield','dilute',6,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

%% boxx 1
zoomb=[-0.35 -0.35 0.7];
bmap='*RdYlBu';
cc= brewermap(8,'Set1'); %distinguishable_colors(3);
%cc=brewermap(3,'Set1');
%cc=distinguishable_color


thick=0.025;

%arrows
arrow(1).start=[-0.26 0.3];
arrow(1).stop=[-0.115 0.25];
arrow(1).faceColor=cc(1,:);

arrow(2).start=[-0.27 0.115];
arrow(2).stop=[-0.15 0.18];
arrow(2).faceColor=cc(2,:); 

arrow(3).start=[-0.24 0.33];
arrow(3).stop=[-0.07 0.3];
arrow(3).faceColor=cc(3,:);

arrow(4).start=[-0.23 0.05];
arrow(4).stop=[-0.125 0.13];
arrow(4).faceColor=cc(4,:);


%drawline
zline=-0.12; 
linStruct.type='--';
linStruct.value=zline;
linStruct.dir='vertical';
% ticks
tic=-0.3:0.1:0.3;

% contours
ctype='dm';


mkmap(1,'type','f',prj,'clims',[-5 5],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','s',prj,'clims',[1.5 3],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','p',prj,'clims',[-1 1.5],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','t',prj,'clims',[7 8],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'drawline',linStruct,'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','vel',prj,'clims',[2.5 3.5],'log','streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'bartag','$\log(v)\,[\mathrm{km/sec}]$','xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','mach',prj,'clims',[0 2.5],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','ztot',prj,'clims',[-1.5 -0.5],'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'contour','dm','contlog','contthick',0.05,'contcolor','k','contlevels',[12:0.5:16],...
'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');%'streamlines',7,

mkmap(1,'type','rho',prj,'clims',[12.5 14],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','n',prj,'clims',[-3.6 -2.2],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','none','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','xray',prj,'clims',[-4 -1.5],'thick',1,...
    'brewer','Greys','labels','none','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

