load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat');

cl=106;
new_env(cl);

prj='XZ';
ind=6;
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
%thick=0.05;
mkmap(4,'type','vel',prj,'clims',[2.5 3.3],'log','vfield','dilute',6,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'bartag','$\log(v)\,[\mathrm{km/sec}]$','print','printag','papIII');

mkmap(4,'type','Mach',prj,'clims',[0 2],'vfield','dilute',6,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
    'print','printag','papIII');

%% boxx 1

%zoomb=[-0.6 -0.3  1.3];
zoomb=[-1 -1  2];
bmap='*RdYlBu';
cc= brewermap(8,'Set1'); %distinguishable_colors(3);
%cc=brewermap(3,'Set1');

thick=0.05; 

%arrows
arrow(1).start=[-0.553 0.868];
arrow(1).stop=[-0.393 0.77];
arrow(1).faceColor=cc(1,:);

arrow(2).start=[0.075 0.855];
arrow(2).stop=[-0.1 0.75];
arrow(2).faceColor=cc(2,:); 

arrow(3).start=[-0.461 0.415];
arrow(3).stop=[-0.326 0.47];
arrow(3).faceColor=cc(3,:);

arrow(4).start=[-0.067 0.472];
arrow(4).stop=[-0.233 0.485];
arrow(4).faceColor=cc(4,:);

arrow(5).start=[-0.31 0.2525];
arrow(5).stop=[-0.2622 0.401];
arrow(5).faceColor=cc(5,:);


%arrow=[];
%drawline
zline=[0.5 0.75]; 
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


mkmap(2,'type','f',prj,'clims',[-5 5],'streamlines',7,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'print','printag','papIII');

mkmap(2,'type','s',prj,'clims',[2 3.5],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'print','printag','papIII');

mkmap(2,'type','p',prj,'clims',[-2 1],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'print','printag','papIII');

mkmap(2,'type','t',prj,'clims',[7 8],'streamlines',7,'thick',thick,...
    'marks',circ(2),'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'drawline',linStruct,'xticks',tic,'yticks',tic,'print','printag','papIII');

mkmap(2,'type','vel',prj,'clims',[2.5 3.5],'log','streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'bartag','$\log(v)\,[\mathrm{km/sec}]$','xticks',tic,'yticks',tic,'print','printag','papIII');

mkmap(2,'type','mach',prj,'clims',[0 2.5],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'print','printag','papIII');

mkmap(2,'type','ztot',prj,'clims',[-2 -0.5],'thick',thick,...
    'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'contour','dm','contlog','contthick',0.05,'contcolor','k','contlevels',[12:0.5:16],...
'xticks',tic,'yticks',tic,'print','printag','papIII');%'streamlines',7,

mkmap(2,'type','rho',prj,'clims',[12.5 14],'streamlines',7,'thick',thick,...
    'brewer',bmap,'labels','full','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'print','printag','papIII');

mkmap(2,'type','xray',prj,'clims',[-4 -2],'streamlines',7,'thick',1,...
    'brewer','Greys','labels','full','title','no','zoom',zoomb,...
    'xticks',tic,'yticks',tic,'print','printag','papIII');

