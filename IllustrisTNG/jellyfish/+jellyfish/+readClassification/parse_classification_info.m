function res = parse_classification_info(cls)
%PARS_CLASSIFICATION parse classification (annotation field in
%classification table) to output the question and answer 
%   Detailed explanation goes here

clsStr=cellstr(cls);
clsStr=string(cls);

res.question=clsStr.extractBetween('"task_label":"','","value"');
res.answer=clsStr.extractBetween('"value":"','"}]');


end
pcolo
