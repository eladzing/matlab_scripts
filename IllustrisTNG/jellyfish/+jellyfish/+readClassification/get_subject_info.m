function res = get_subject_info(subject_data)
%GET_SUBJECT_INFO parse the subject_data field form classification table
%and return the important info
%   given the subject data field from the classification table, which is a
%   STRING class object, parse the srting and extract relevant information:
%   simulation, snapshot, host ID and Subfind ID

res.subjectID=str2double(subject_data.extractBetween('{"','":{"ret'));

simname0=subject_data.extractBetween('"Filename":"',"_");

%nval=simname0.extractBetween('n','TNG');
switch simname0.extractBetween('L','n')
    case('205')
        switch simname0.extractBetween('n','TNG')
            case '2500'
                res.simname='TNG300';
            case '1250'
                res.simname='TNG300-2';
            case '625'
                res.simname='TNG300-3';
            otherwise
                error('get_subjectinfo - unknown resolution: %s', simname0)
        end
    case('75')
        switch simname0.extractBetween('n','TNG')
            case '1820'
                res.simname='TNG100';
            case '910'
                res.simname='TNG100-2';
            case '455'
                res.simname='TNG100-3';
            otherwise
                error('get_subjectinfo - unknown resolution: %s', simname0)
        end
        res.simname='TNG100';
    case('35')
        switch simname0.extractBetween('n','TNG')
            case '2160'
                res.simname='TNG50';
            case '1080'
                res.simname='TNG50-2';
            case '540'
                res.simname='TNG50-3';
            case '270'
                res.simname='TNG50-4';
            otherwise
                error('get_subjectinfo - unknown resolution: %s', simname0)
        end
        
    otherwise
        res.simname='nonTNG';
        %error('get_subjectinfo - Illegal simulation name: %s', simname0)
end

res.fullSimName=simname0;


snap0=subject_data.extractBetween('TNG',"_HostFofID");
res.snap=str2double(snap0.extractAfter('_0'));

res.hostID=str2double(subject_data.extractBetween("HostFofID_","_SubfindID"));

res.subfindID=str2double(subject_data.extractBetween("SubfindID_",".png"));


end

