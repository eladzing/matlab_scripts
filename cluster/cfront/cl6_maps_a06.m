load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat');

cl=6;
new_env(cl,'a06');
cc= brewermap(8,'Set1');

prj='YZ';
ind=10;
bmap='*RdYlBu';

arrow(1).start=[-0.17 0.277];
arrow(1).stop=[-0.07 0.221 ];
arrow(1).faceColor=cc(1,:);

arrow(2).start=[-0.15 -0.383];
arrow(2).stop=[-0.037 -0.35];
arrow(2).faceColor=cc(5,:);

arrow(3).start=[-0.127 0.107];
arrow(3).stop=[-0.058 0.157];
arrow(3).faceColor=cc(2,:);

arrow(4).start=[-0.192 -0.192];
arrow(4).stop=[-0.106 -0.233 ];
arrow(4).faceColor=cc(3,:);

arrow(5).start=[-0.143 -0.105];
arrow(5).stop=[-0.028 -0.14];
arrow(5).faceColor=cc(4,:);

zoomb=[-0.4 -0.4 0.8];thick=0.2;
mkmap(2,'type','t',prj,'clims',[7 8],'streamlines',7,'proper','thick',thick,...
'brewer',bmap,'labels','full','title','no',...
'zoom',zoomb,'print','printag','papIII','arrow',arrow);

zoomb=[-0.4 -0.4 0.8];thick=0.2;
mkmap(2,'type','gas',prj,'clims',[13 14.5],'streamlines',7,'proper','thick',thick,...
'brewer',bmap,'labels','full','title','no',...
'zoom',zoomb,'print','printag','papIII','arrow',arrow);

zoomb=[-0.4 -0.4 0.8];thick=0.2;
mkmap(2,'type','n',prj,'clims',[-3.2 -1.6],'streamlines',7,'proper','thick',thick,...
'brewer',bmap,'labels','full','title','no',...
'zoom',zoomb,'print','printag','papIII','arrow',arrow);

thick=0.2;
zoomb=[-1 -1 2];
mkmap(4,'type','f',prj,'clims',[-3 3],'vfield','dilute',8,'proper','thick',thick,...
'brewer',bmap,'labels','full','title','no',...
'print','printag','papIII','zoom',zoomb);


% shk=shockedge_a06;
% shkM1=shockedgeMask1_a06;
% shkM2=shockedgeMask2_a06;
% c1=shk{ind,4};
% c2=shk{ind,5};
% circ=cat(2,c1(shkM1(ind,:)),c2(shkM2(ind,:)));

% % boxx 8
% thick=0.5;
% 
% mkmap(8,'type','rho',prj,'clims',[7 14],'vfield','dilute',6,'proper','thick',thick,...
%     'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
%     'print','printag','papIII');
% 
% mkmap(8,'type','dm',prj,'clims',[7 15],'proper','thick',thick,...
%     'marks',circ(2),'brewer',bmap,'labels','full','title','no',...
%     'print','printag','papIII');

%%boxx 4 
