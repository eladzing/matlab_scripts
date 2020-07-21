function res = simulationMask(clsTab,simName)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


clsTab(~ismember(clsTab.simName,simName),:)=[];

res=clsTab;

end

