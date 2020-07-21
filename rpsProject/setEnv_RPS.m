function  setEnv_RPS()
%SETENV_RPS set the environment for RPS project 

global DEFAULT_PRINTOUT_DIR
global DEFAULT_FIG_DIR
global DEFAULT_MOVIE_DIR
global DEFAULT_MATFILE_DIR
global DEFAULT_TNG_MATFILE_DIR
global DRACOFLAG


HOME_PATH='/home/zinger/workProjects/clusterRPS';
DRACO_PATH='/isaac/ptmp/gc/eladzing/RPS_project';

DRACO_PRINTOUT_DIR = [DRACO_PATH '/printout'];
HOME_PRINTOUT_DIR = [HOME_PATH '/printout'];

DRACO_MATFILE_DIR = [DRACO_PATH '/matFiles'];
HOME_MATFILE_DIR = '/home/zinger/workProjects/matlab_scripts/rpsProject/matFiles';

HOME_TNG_MATFILE_DIR  = '/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles';
DRACO_TNG_MATFILE_DIR = '/isaac/ptmp/gc/eladzing/IllustrisTNG/matFiles';



%% Check to see if we are in draco

pp=pwd;
DRACOFLAG=strcmp(pp(1:7),'/draco/')  || strcmp(pp(1:7),'/isaac/');


if ~DRACOFLAG
% set paths


DEFAULT_PRINTOUT_DIR=HOME_PRINTOUT_DIR;
DEFAULT_MATFILE_DIR=HOME_MATFILE_DIR;
DEFAULT_TNG_MATFILE_DIR=HOME_TNG_MATFILE_DIR ;

else
    
DEFAULT_PRINTOUT_DIR=DRACO_PRINTOUT_DIR;
DEFAULT_MATFILE_DIR=DRACO_MATFILE_DIR;
DEFAULT_TNG_MATFILE_DIR =DRACO_TNG_MATFILE_DIR ;
end

DEFAULT_FIG_DIR=[DEFAULT_PRINTOUT_DIR '/figFiles'];
DEFAULT_MOVIE_DIR=[DEFAULT_PRINTOUT_DIR '/movies'];


global cosmoStruct

cosmoStruct=set_LCDM_cosmology;

end

