function cluster(clst,a)

%global FILE_FORMAT;
%base_dir = '/home/alf/eladzing/data/kravtsov/clusters';
%inpath= sprintf('%s/%s',base_dir,clst);
%filform= '%s_a%sL%dMpc.dat';
%filform=sprintf(filform,clst,a);

FILE_FORMAT = sprintf(FILE_FORMAT,clst,clst,a);
