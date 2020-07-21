function fofStruct =  addTvirFofs( fofStruct)
%ADDTVIRFOFS -  calculate TVIR for the FOFs and add as fields to structure. 

global illUnits
global cosmoStruct
% M_Mean200
[~, ~, tvir, ~]=calculate_virials('mvir',fofStruct.Group_M_Mean200*illUnits.massUnit,'cosmo',cosmoStruct,'delv',200,'mean');
fofStruct.Group_T_Mean200=tvir;

% M_Crit200
[~, ~, tvir, ~]=calculate_virials('mvir',fofStruct.Group_M_Crit200*illUnits.massUnit,'cosmo',cosmoStruct,'delv',200,'crit');
fofStruct.Group_T_Crit200=tvir;

% M_Crit500
[~, ~, tvir, ~]=calculate_virials('mvir',fofStruct.Group_M_Crit500*illUnits.massUnit,'cosmo',cosmoStruct,'delv',500,'crit');
fofStruct.Group_T_Crit500=tvir;

% M_TopHat200
[~, ~, tvir, ~]=calculate_virials('mvir',fofStruct.Group_M_TopHat200*illUnits.massUnit,'cosmo',cosmoStruct,'crit');
fofStruct.Group_T_TopHat200=tvir;


end

