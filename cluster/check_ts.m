list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    
tmins=[];
bigbox=8;
global FILE_FORMAT_SPHERE;
for id=1:length(list1)
    halopath=sprintf(hpath,list1(id));
    clustername=sprintf('CL%d',list1(id))
    switch list1(id)
        case {101,102,103,104,105,106,107,5}
            smallbox=2;
        otherwise
            smallbox=1;
    end
    FILE_FORMAT_SPHERE = sprintf('%s/%s', halopath, '%s_sphere_%d.mat');
    %for ii=(log2(smallbox)):log2(bigbox)
        boxx=smallbox
        ts =  T_sphere(boxx);
        ts = ts(1:255,:,:);
        tmins(id,ii+1)=min(ts(:));
    %end
     
end

