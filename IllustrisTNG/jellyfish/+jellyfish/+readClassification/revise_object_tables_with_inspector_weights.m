
pilotTag(length(newScore.ncls0),1)="Why";
pilotTag(newScore.ncls0==0)="No";
pilotTag(newScore.ncls0==5)="TNG100";
pilotTag(newScore.ncls0==6)="TNG50";





objectTable=table(objectCJF.subject_ids,objectCJF.sim,objectCJF.snap,...
    objectCJF.subfind,objectCJF.type,objectCJF.tag,...
    objectCJF.score,objectCJF.scoreTotal,objectCJF.clsNum,...
    newScore.score_hi_down_up_exp',newScore.weight',newScore.ncls',...
    (newScore.nexp'+newScore.ncls0'),pilotTag,...
    'variableNames',{'subject_ids','sim','snap','subfind','type','tag',...
    'scoreRaw','scoreTotalRaw','clsNumRaw',...
    'scoreWeighted','weight','clsNumWeighted','expNum','pilot'});
msk=objectTable.clsNumWeighted==0;
objectTable.scoreWeighted(msk)=0;

objectTablePref=table(objectTablePrefOriginal.subject_ids,objectTablePrefOriginal.sim,objectTablePrefOriginal.snap,...
    objectTablePrefOriginal.subfind,objectTablePrefOriginal.type,objectTablePrefOriginal.tag,...
    objectTablePrefOriginal.score,objectTablePrefOriginal.scoreTotal,objectTablePrefOriginal.clsNum,...
    newScorePref.score_hi_down_up_exp',newScorePref.weight',newScorePref.ncls',...
    newScorePref.nexp',...
    'variableNames',{'subject_ids','sim','snap','subfind','type','tag',...
    'scoreRaw','scoreTotalRaw','clsNumRaw',...
    'scoreWeighted','weight','clsNumWeighted','expNum'});

len=height(objectTablePref);
scoreRand=zeros(len,1);
weightRand=zeros(len,1);
for i=1:len
    searchTag=objectTablePref.tag(i).replace("pref","rand");  
    ind=find(objectTable.tag==searchTag);
    if isempty(ind)
        error('no match %i',i)
    end
    scoreRand(i)=objectTable.scoreWeighted(ind);
    weightRand(i)=objectTable.weight(ind);
end


objectTableComp=table(objectTableCompOriginal.subject_ids,objectTableCompOriginal.sim,...
    objectTableCompOriginal.snap,objectTableCompOriginal.subfind,objectTableCompOriginal.tag,...
    objectTableCompOriginal.scorePref,objectTableCompOriginal.scoreRand,...
     newScorePref.score_hi_down_up_exp',scoreRand,...
     newScorePref.weight',weightRand,...
     'variableNames',{'subject_ids','sim','snap','subfind','tag',...
    'scoreRawPref','scoreRawRand','scoreWeightedPref','scoreWeightedRand',...
    'weightPref','weightRand'});



%% 
cnt=0;
for i=1:length(indx20)
    searchTag=objectTable.tag(indx20(i)).replace("rand","pref");  
    ind=find(objectTableComp.tag==searchTag);
    if ~isempty(ind)
       cnt=cnt+1;
    end
    
end