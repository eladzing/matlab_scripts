function res = isJellyfish_single(cls,checkQuestion)
%ISJELLYFISH_SINGLE returns a logical true/false for a single
%calssification event (recorded in 'annotation' field
%   The input data is of type CATEGORICAL' and is parsed into a question
%   and answer. First check that its the right question and then check the
%   answer.

checkFlag=false;

if exist('checkQuestion','var')
    checkFlag=strcmp(checkQuestion,'check');
end

rightQuestion='Do you think that the galaxy at the center looks like a jellyfish ?';


q_and_a=parse_classification_info(cls);

if checkFlag && ~strcmp(q_and_a.question,rightQuestion)
        error('Wrong classification question: %s',q_and_a.question);
end

res=strcmp(q_and_a.answer,"Yes");


end

