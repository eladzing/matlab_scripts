
a='a06';
find_r500c;

cl=9;
new_env(cl,a)
ii=find(cl==clust);

%r5=r500c(ii).*get_rvir;

r5=r500c(ii).*hub;

zmc=ceil(100*r5*1.05)/100;
boxx=2.^ceil(log(2*r5*(1+zred))./log(2));
zmc=min(zmc,boxx./(1+zred)/2);
%zmc=boxx/(1+zred)/2;
zm=[-1 -1 2].*zmc;

% mkmap(boxx,'type','xray','vfield','dilute',8,'contlevels',-4.5:0.5:0,...
%     'contourtype','xray','contthick',boxx,'marks',r500c(ii),'markcolor','w',...
%     'zoom',zm,'thick',boxx,'contcolor','k','grid','proper')

mkmap(boxx,'type','f','clims',[-3 3],'streamlines',3,'thick',0.1,...
    'contourtype','xray','contthick',boxx,'contcolor','k','contlevels',-4.5:0.5:0,...
    'marks',r500c(ii),'markcolor','w','zoom',zm,'proper',...
    'title','none','print','printag','mrgTest','labels','half')


tot=RHOTOT(boxx);
tot=smooth3(tot,'box',3);
mkmap(boxx,'data',tot,'vfield','dilute',4,'clims',[12.5 15],'log','thick',boxx,...
    'contourtype','xray','contthick',boxx,'contcolor','k','contlevels',6:0.5:12.5,...
    'marks',[r500c(ii) get_rvir*0.25],'zoom',zm,'proper','grid')