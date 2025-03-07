function res=format_classification_table(clsTab)


% focus on the right workflow
right_wf=11990;
clsTab(~ismember(clsTab.workflow_id,right_wf),:)=[];

% remove unwanted columns 
clsTab=removevars(clsTab,{'user_ip','workflow_name','workflow_version',...
    'metadata','gold_standard','expert','workflow_id'});
    
% add snapshot as column 
clsTab.snap=clsTab.subject_data.extractBetween('TNG_',"_HostFofID");
clsTab.snap=clsTab.snap.extractAfter('0').double;

simBase=clsTab.subject_data.extractBetween('"Filename":"',"_");
sim=simBase.extractBetween('L','n');
resol=simBase.extractBetween('n','TNG');

sim(strcmp(sim,'35') & strcmp(resol,'2160'))='TNG50';
sim(strcmp(sim,'35') & strcmp(resol,'1080'))='TNG50-2';
sim(strcmp(sim,'35') & strcmp(resol,'540'))='TNG50-3';
sim(strcmp(sim,'35') & strcmp(resol,'270'))='TNG50-4';

sim(strcmp(sim,'75') & strcmp(resol,'1820'))='TNG100';
sim(strcmp(sim,'75') & strcmp(resol,'910'))='TNG100-2';
sim(strcmp(sim,'75') & strcmp(resol,'455'))='TNG100-3';

sim(strcmp(sim,'205') & strcmp(resol,'2500'))='TNG300';
sim(strcmp(sim,'205') & strcmp(resol,'1250'))='TNG300-2';
sim(strcmp(sim,'205') & strcmp(resol,'625'))='TNG300-3';

clsTab.simName=sim;

clsTab.host=clsTab.subject_data.extractBetween("HostFofID_","_SubfindID").double;
clsTab.subfind=clsTab.subject_data.extractBetween("SubfindID_",".png").double;

% replace hyphens with underscores in user_names
clsTab.user_name=clsTab.user_name.replace("-","_");

%% get answer 
%clsTab.question=string(clsTab.annotations).extractBetween('"task_label":"','","value"');
clsTab.isJelly=strcmp(string(clsTab.annotations).extractBetween('"value":"','"}]'),"Yes");

%% cleanup
clsTab=removevars(clsTab,{'subject_data','annotations'});

res=clsTab;

