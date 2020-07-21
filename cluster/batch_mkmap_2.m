% % cl6 
 new_env(6,'csf','a1');

% flux 
% mkmap(8,'marks',get_rvir.*[0.01 0.2],'type','density','proj',[1 0 0],'vfield','dilute',6,'title','off',...
%     'labels','units','print','format','both','printag','struct')
% mkmap(8,'marks',get_rvir.*[0.01 0.2],'type','entropy','proj',[1 0 0],'vfield','dilute',6,'title','off',...
%     'labels','units','print','format','both','printag','struct','clims',[-1 2])
% mkmap(8,'marks',get_rvir.*[0.01 0.2],'type','flux','proj',[1 0 0],'vfield','dilute',6,'title','off',...
%     'labels','units','print','format','both','printag','struct','clims',[-3 3])
mkmap(4,'marks',get_rvir.*[0.01 0.2],'type','flux','proj',[1 0 0],'vfield','dilute',6,'title','off',...
    'labels','units','print','format','both','printag','struct','clims',[-3 3])
mkmap(2,'marks',get_rvir.*[0.01 0.2],'type','flux','proj',[1 0 0],'vfield','dilute',6,'title','off',...
    'labels','units','print','format','both','printag','struct','clims',[-3 3])
% mkmap(8,'marks',get_rvir.*[0.01 0.2],'type','temperature','proj',[1 0 0],'vfield','dilute',6,'title','off',...
%     'labels','units','print','format','both','printag','struct','clims',[4.5 8])
% 
% % % cl3
% new_env(3,'csf','a1');
% 
% mkmap(4,'marks',get_rvir.*[0.01 0.2],'type','flux','proj',[0 0 1],'vfield','dilute',6,'title','off',...
%     'labels','units','print','format','both','printag','struct','clims',[-3 3])
% mkmap(2,'marks',get_rvir.*[0.01 0.2],'type','flux','proj',[0 0 1],'vfield','dilute',6,'title','off',...
%     'labels','units','print','format','both','printag','struct','clims',[-3 3])
% 
% % cl5
% new_env(5,'csf','a1');
% global zred
% mkmap(4,'marks',get_rvir.*[0.01 0.2],'type','flux','proj','xz','vfield','dilute',6,'titletype','custom','title',sprintf('$z=%3.2g$',zred),...
%     'labels','units','print','format','both','printag','struct','clims',[-3 3])
% 
% new_env(5,'csf','a06');
% mkmap(4,'marks',get_rvir.*[0.01 0.2],'type','flux','proj','xz','vfield','dilute',6,'titletype','custom','title',sprintf('$z=%3.2g$',zred),...
%     'labels','units','print','format','both','printag','struct','clims',[-3 3],'proper')
% 
% % cl11
% new_env(11,'csf','a1');
% mkmap(2,'marks',get_rvir.*[0.01 0.2],'type','flux','proj','xy','vfield','dilute',6,'titletype','custom','title',sprintf('$z=%3.2g$',zred),...
%     'labels','units','print','format','both','printag','struct','clims',[-3 3])
% 
% new_env(11,'csf','a06');
% mkmap(2,'marks',get_rvir.*[0.01 0.2],'type','flux','proj','xy','vfield','dilute',6,'titletype','custom','title',sprintf('$z=%3.2g$',zred),...
%     'labels','units','print','format','both','printag','struct','clims',[-3 3],'proper')
