% read in table and format it to our liking.
clsTab=jellyfish.readClassification.import_classification_file_CJF('191121');

%% separate the Phase 1 classification table
startDate=datetime('01/06/21 00:00','InputFormat',"dd/MM/yy HH:mm");
endDate=datetime('14/06/21 23:59','InputFormat',"dd/MM/yy HH:mm");

clsTabP1=jellyfish.readClassification.format_classification_table_CJF(clsTab,'start',startDate,'end',endDate);
%head(clsTab,5)
%clsTab=jellyfish.readClassification.simulationMask(clsTab,'TNG50');
% Classification info based on data downloaded on 14/7/21
% review and training set objecs classifications are removed. 
% The raw classification data is reformatted into a new table with columns contiaing the snapshot, id  Id for the object as well as the label (jelly or not):

%% Separaete the phase 2 classification table

startDate=datetime('22/08/21 00:00','InputFormat',"dd/MM/yy HH:mm");
endDate=datetime('19/11/21 23:59','InputFormat',"dd/MM/yy HH:mm");

clsTabP2=jellyfish.readClassification.format_classification_table_CJF(clsTab,'start',startDate,'end',endDate);

%% Generate object tables for the two phases 

% number of objects
[idListP1, ~,~]=unique(clsTabP1.subject_ids);
ngalP1=length(idListP1);

fprintf('there are %i objects in Phase I\n',ngalP1)

[idListP2,~,~]=unique(clsTabP2.subject_ids);
ngalP2=length(idListP2);
fprintf('there are %i objects in Phase II\n',ngalP2)


% create object table
objectTableP1=jellyfish.readClassification.build_object_table(clsTabP1);
objectTableP2=jellyfish.readClassification.build_object_table(clsTabP2);
objectTableP2=objectTableP2(5:end,:);  % remove first 4 objects with only 1 classification 
% table(idList,clsTab.simName(ia),clsTab.snap(ia),clsTab.subfind(ia),zeros(size(idList)),zeros(size(idList)),...
%     'variableNames',{'subject_ids','sim','snap','subfind','Ncls','Ycls'});

ngalP1=height(objectTableP1);
ngalP2=height(objectTableP2);
fprintf('there are %i objects in Phase I\n',ngalP1)
fprintf('there are %i objects in Phase II\n',ngalP2)


snapsP1=unique(objectTableP1.snap);
zredsP1=round(10.*illustris.utils.snap2redshift(snapsP1))./10;

snapsP2=unique(objectTableP2.snap);
zredsP2=round(10.*illustris.utils.snap2redshift(snapsP2))./10;

%% Identify doubles objects in phase II 

doubleIndP1=[];
doubleIndP2=[];

cnt=0;
for i=1:height(objectTableP1)
    
%     if mod(i,1000)==0
%         fprintf('%i ',i/1000);
%     end
    
    msk=strcmp(objectTableP1.sim(i),objectTableP2.sim) & objectTableP2.snap==objectTableP1.snap(i);
    ii=find(msk & objectTableP2.subfind==objectTableP1.subfind(i));
    
    if ~isempty(ii)
        cnt=cnt+1;
        doubleIndP1(cnt)=i;  % index of object in table 
        doubleIndP2(cnt)=ii; % index of object in table 
        
    end
    
end

dp1=objectTableP1.score(doubleIndP1);
dp2=objectTableP2.score(doubleIndP2);

%% stitch togehter both object table after removing doubles 
msk=true(height(objectTableP2),1);
msk(doubleIndP2)=false;
tempTab=objectTableP2(msk,:); % nondouble copy of the phase II table 

objectTable=[objectTableP1 ; tempTab];

%% generate double tables 
doubleTable=objectTableP1(doubleIndP1,1:6);

doubleTable.Properties.VariableNames(6)={'scoreI'};
doubleTable.scoreII=objectTableP2.score(doubleIndP2);

%% save object table 






