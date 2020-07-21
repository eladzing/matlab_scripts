list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat');

cl=24;
prj='xy';

ind=find(cl==list,1,'first');

new_env(cl);

dmap='*RdBu';
smap='YlOrRd';


% shk=shockedge_a1;
% shkM1=shockedgeMask1_a1;
% shkM2=shockedgeMask2_a1;
% c1=shk{ind,8};
% c2=shk{ind,5};
% circ=cat(2,c1(shkM1(ind,:)),c2(shkM2(ind,:)));

% boxx 8
thick=0.5;
mkmap(8,'type','f',prj,'clims',[-3 3],'vfield','dilute',6,'thick',thick,...
    'brewer',dmap,'labels','full','title','no');
%
% mkmap(8,'type','rho',prj,'clims',[7 14],'vfield','dilute',6,'thick',thick,...
%     'marks',circ(2),'brewer',bmap,'labels','full','title','no');
%
% mkmap(8,'type','dm',prj,'clims',[7 15],'thick',thick,...
%     'marks',circ(2),'brewer',bmap,'labels','full','title','no');
%
%%boxx 4
thick=0.100;
mkmap(4,'type','f',prj,'clims',[-3 3],'vfield','dilute',6,'thick',thick,...
    'brewer',dmap,'labels','full','title','no');
mkmap(4,'type','s',prj,'clims',[1.5 3],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no');

mkmap(4,'type','p',prj,'clims',[-1 1.5],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no');

mkmap(4,'type','t',prj,'clims',[7 8],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no');

mkmap(4,'type','rho',prj,'clims',[12.5 14],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no');
%%boxx 2
%bmap='*Spectral';
thick=0.05;
mkmap(2,'type','vel','norm',get_vvir,prj,'clims',[-1 0.5],'log','vfield','dilute',6,'thick',thick,...
    'brewer',smap,'labels','full','title','no',...
    'bartag','$\log(v/V_{\mathrm{vir}})$');

mkmap(2,'type','Mach',prj,'clims',[0 2],'vfield','dilute',6,'thick',thick,...
    'brewer',smap,'labels','full','title','no');

mkmap(2,'type','flux',prj,'clims',[-3 3],'vfield','dilute',6,'thick',thick,...
    'brewer',dmap,'labels','full','title','no','marks',0.25*get_rvir);


mkmap(2,'type','s',prj,'clims',[1.5 3],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no');

mkmap(2,'type','p',prj,'clims',[-1 1.5],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no');

mkmap(2,'type','t',prj,'clims',[7 8],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no');

mkmap(2,'type','rho',prj,'clims',[12.5 14],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no');

mkmap(2,'type','ztot',prj,'clims',[-1.5 -0.5],'thick',thick,...
    'brewer',smap,'labels','full','title','no',...
    'contour','dm','contlog','contthick',0.05,'contcolor','k','contlevels',12:0.5:16);%'streamlines',7,

%% boxx 1
zoomb=[-0.5 -0.5 1.0];
%zoomb=[-0.35 -0.35 0.7];
bmap='*RdYlBu';
cc= brewermap(8,'Set1'); %distinguishable_colors(3);
%cc=brewermap(3,'Set1');
%cc=distinguishable_color


thick=0.025;


% contours
ctype='dm';


mkmap(1,'type','f',prj,'clims',[-5 5],'streamlines',7,'thick',thick,...
    'brewer',dmap,'labels','full','title','no','zoom',zoomb,...
    'marks',0.25*get_rvir);

mkmap(1,'type','s',prj,'clims',[1.5 3],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no','zoom',zoomb);

mkmap(1,'type','p',prj,'clims',[-1 1.5],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no','zoom',zoomb);

mkmap(1,'type','t',prj,'clims',[7 8],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no','zoom',zoomb);

% mkmap(1,'type','vel',prj,'clims',[2.5 3.5],'log','streamlines',7,'thick',thick,...
%     'brewer',smap,'labels','full','title','no','zoom',zoomb,...
%     'bartag','$\log(v)\,[\mathrm{km/sec}]$','xticks',tic,'yticks',tic,'marks',0.25*get_rvir);
% 
% mkmap(1,'type','mach',prj,'clims',[0 2.5],'streamlines',7,'thick',thick,...
%     'brewer',smap,'labels','full','title','no','zoom',zoomb,...
%     'xticks',tic,'yticks',tic);

mkmap(1,'type','ztot',prj,'clims',[-1.5 -0.5],'thick',thick,...
    'brewer',smap,'labels','full','title','no','zoom',zoomb,...
    'contour','dm','contlog','contthick',0.05,'contcolor','k','contlevels',12:0.5:16)
    

mkmap(1,'type','rho',prj,'clims',[12.5 14],'streamlines',7,'thick',thick,...
    'brewer',smap,'labels','full','title','no','zoom',zoomb);

% mkmap(1,'type','xray',prj,'clims',[-4 -1.5],'thick',1,...
%     'brewer','Greys','labels','full','title','no','zoom',zoomb,...
%     'xticks',tic,'yticks',tic);
% 
