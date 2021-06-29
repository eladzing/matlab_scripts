function cluster(clst,a)
env;
global FILE_FORMAT;
global CLUSTER_PATH
BASE_DIR='/home/alf/eladzing/data/kravtsov/clusters';
CLUSTER_PATH=sprintf('%s/%s',BASE_DIR,clst);
inpath= sprintf('%s/%s',CLUSTER_PATH,FILE_FORMAT);
%filform= '%s_a%sL%dMpc.dat';
%filform=sprintf(filform,clst,a);

FILE_FORMAT = sprintf(inpath,'%s',a,'%d');
