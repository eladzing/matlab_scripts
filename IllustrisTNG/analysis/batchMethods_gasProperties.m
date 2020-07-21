
%% run over methods boxes and create gas  properties 

%method_list={'0000','2201','3000','3101','3102','3103','3104',...
%    '3301','3302','3403','3404','3801','3802','2101','2302','0030'};

method_list={'2302','0030'};


snap=4;
for i=1:length(method_list)
    
    bp=illustris.set_env(35,method_list{i});
    readFlag=true;
    
    generate_gasProperties_fofs;
    generate_gasProperties_light;
    generate_bh_gals
    
    
    fprintf('=========================================================================== \n');
end