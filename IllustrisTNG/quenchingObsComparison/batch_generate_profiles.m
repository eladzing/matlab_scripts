%% batch generate profiles 
list={'100','100-2','100-3',...
    '300','300-2','300-3',...
    '50','50-2','50-3','50-4'};

snap=99;
readFlag=true;
for i=1:length(list)
    
    if contains(list{i},'50')
        massThresh=10^7;
    else
        massThresh=10^9;
    end
    
    bp=illustris.set_env(list{i});
    
    %comparison_take_2
    projected_profiles
    
    clear profStruct
end
