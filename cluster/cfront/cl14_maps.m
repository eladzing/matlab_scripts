load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat');

cl=14;
new_env(cl,'a06');

prj='XZ';
ind=15;
bmap='*RdYlBu';


shk=shockedge_a06;
shkM1=shockedgeMask1_a06;
shkM2=shockedgeMask2_a06;
c1=shk{ind,4};
c2=shk{ind,5};
circ=cat(2,c1(shkM1(ind,:)),c2(shkM2(ind,:)));

% boxx 8
thick=0.05;

mkmap(2,'type','rho',prj,'clims',[7 14],'vfield','dilute',6,'thick',thick,...
    'marks',circ(3),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

mkmap(2,'type','dm',prj,'clims',[7 15],'thick',thick,...
    'marks',circ(3),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

%%boxx 4 
%thick=0.100;
mkmap(2,'type','f',prj,'clims',[-3 3],'vfield','dilute',6,'thick',thick,...
    'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

%%boxx 2
bmap='*Spectral';
%thick=0.05;
mkmap(2,'type','vel',prj,'clims',[1.8 3.3],'log','vfield','dilute',6,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

mkmap(2,'type','Mach',prj,'clims',[0 2],'vfield','dilute',6,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

%% boxx 1
zoomb=[-0.24 -0.18  0.48];
bmap='*RdYlBu';
cc= brewermap(8,'Set1'); %distinguishable_colors(3);
%cc=brewermap(3,'Set1');

thick=0.025; 

%arrows
arrow(1).start=[-0.553 0.868];
arrow(1).stop=[-0.393 0.77];
arrow(1).faceColor='r';

arrow(2).start=[0.075 0.855];
arrow(2).stop=[-0.1 0.75];
arrow(2).faceColor='b'; 

arrow(3).start=[-0.461 0.415];
arrow(3).stop=[-0.326 0.47];
arrow(3).faceColor='r';

arrow(4).start=[-0.067 0.446];
arrow(4).stop=[-0.233 0.485];
arrow(4).faceColor='b';


%drawline
zline=[-0.25 -0.120 ]; 
linStruct(1).type='--';
linStruct(1).value=zline(1);
linStruct(1).dir='horizontal';
linStruct(2).type='--';
linStruct(2).value=zline(2);
linStruct(2).dir='horizontal';
% ticks
tic=-0.6:0.2:1;


% contours
ctype='dm';


mkmap(1,'type','f',prj,'clims',[-5 5],'streamlines',7,'thick',thick,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,'proper');
    %'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','s',prj,'clims',[2 3.5],'streamlines',7,'thick',thick,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,'proper');
    %'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','p',prj,'clims',[-2 1],'streamlines',7,'thick',thick,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,'proper');
    %'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','t',prj,'clims',[7 8],'streamlines',7,'thick',thick,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,'proper');
    %'drawline',linStruct,'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','vel',prj,'clims',[2.5 3.5],'log','streamlines',7,'thick',thick,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,'proper');
    %'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','mach',prj,'clims',[0 2.5],'streamlines',7,'thick',thick,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,'proper');
    %'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','ztot',prj,'clims',[-2 -0.5],'thick',thick,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'contour','dm','contlog','contthick',0.05,'contcolor','k','contlevels',12:0.5:16,'proper');
%'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');%'streamlines',7,

mkmap(1,'type','rho',prj,'clims',[12.5 14],'streamlines',7,'thick',thick,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,'proper');
 %   'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

mkmap(1,'type','xray',prj,'clims',[-4 -2],'streamlines',7,'thick',1,...
   'brewer',bmap,'labels','full','title','no','zoom',zoomb,'proper');
  %  'xticks',tic,'yticks',tic,'arrow',arrow,'print','printag','papIII');

