function res=format_classification_table_CJF(clsTab,varargin)
%% format the classification table
% option for trimming by dates. 

%% parse arguments 

i=1;
startDate=datetime('01/01/00 00:00','InputFormat',"dd/MM/yy HH:mm");
endDate=datetime('31/12/40 23:59','InputFormat',"dd/MM/yy HH:mm");

while(i<=length(varargin))
    switch(lower(varargin{i}))
        case {'start','starts','startdate','from'}
            i=i+1;
            dt=varargin{i};
            if isdatetime(dt)
                startDate=dt;
            else
                error('%s - start date must be a datetime object',current_function);
            end
            
        case {'end','ends','enddate','to','until','till'}
            i=i+1;
            dt=varargin{i};
            if isdatetime(dt)
                endDate=dt;
            else
                error('%s - start date must be a datetime object',current_function);
            end
            
        otherwise
            error('%s - Illegal argument: %s',current_function,varargin{i});
    end
    i=i+1;
end



% focus on the right workflow - not relavent anymore 
% right_wf=11990;
% clsTab(~ismember(clsTab.workflow_id,right_wf),:)=[];

%% remove unwanted columns 
clsTab=removevars(clsTab,{'user_ip','workflow_name','workflow_version',...
    'metadata','gold_standard','expert','workflow_id'});
    
%% remove classification before startdate and after end date 
% replace 'created at' with a datetime. 
clsTab.created_at=datetime(clsTab.created_at,'InputFormat','yyyy-MM-dd HH:mm:ss ''UTC');

dateMask=clsTab.created_at>=startDate & clsTab.created_at<=endDate;

clsTab=clsTab(dateMask,:);


%% add snapshot as column 
clsTab.snap=clsTab.subject_data.extractAfter('_snap');
clsTab.snap=clsTab.snap.extractBefore('_');
clsTab.snap=clsTab.snap.extractAfter('0').double;

simBase=clsTab.subject_data.lower.extractAfter('filename":"');
simBase1=clsTab.subject_data.lower.extractAfter('filename1":"');
simBase(ismissing(simBase))=simBase1(~ismissing(simBase1));

sim=simBase.extractBefore('n');

resol=simBase.extractBefore('tng');
resol=resol.extractAfter('n');




sim(strcmp(sim,'l35') & strcmp(resol,'2160'))='TNG50';
sim(strcmp(sim,'l35') & strcmp(resol,'1080'))='TNG50-2';
sim(strcmp(sim,'l35') & strcmp(resol,'540'))='TNG50-3';
sim(strcmp(sim,'l35') & strcmp(resol,'270'))='TNG50-4';

sim(strcmp(sim,'l75') & strcmp(resol,'1820'))='TNG100';
sim(strcmp(sim,'l75') & strcmp(resol,'910'))='TNG100-2';
sim(strcmp(sim,'l75') & strcmp(resol,'455'))='TNG100-3';

sim(strcmp(sim,'l205') & strcmp(resol,'2500'))='TNG300';
sim(strcmp(sim,'l205') & strcmp(resol,'1250'))='TNG300-2';
sim(strcmp(sim,'l205') & strcmp(resol,'625'))='TNG300-3';

clsTab.simName=sim;

%clsTab.host=clsTab.subject_data.extractBetween("HostFofID_","_SubfindID").double;
subhalo=simBase.extractAfter("subhalo");
clsTab.subfind=subhalo.extractBefore('_snap').double;

% replace hyphens with underscores in user_names
clsTab.user_name=clsTab.user_name.replace("-","_");


%% 
%% get answer 
%clsTab.question=string(clsTab.annotations).extractBetween('"task_label":"','","value"');
answer=string(clsTab.annotations).lower.extractAfter('value');
answer=answer.extractBefore('"}');
clsTab.isJelly=answer.contains('yes');%  strcmp(string(clsTab.annotations).extractBetween('"value":"','"}]'),"Yes");

%% cleanup
clsTab=removevars(clsTab,{'subject_data','annotations'});

%% remove training & review
removalMask=simBase.contains('_training_sample_') | simBase.contains('_review_sample_');
clsTab(removalMask,:)=[];

res=clsTab;

