function basePath=set_env(varargin)

% defaults
global simDisplayName
global simName
global LBox
global DEFAULT_PRINTOUT_DIR
global DEFAULT_MOVIE_DIR
global POSTPROCESSING
global DEFAULT_MATFILE_DIR
global DEFAULT_REMOTE_MATFILE_DIR
%global TEMPFILEDIR
global BASEPATH
global DRACOFLAG
global DEFAULT_FIG_DIR
global myPOSTPROCESSING
global JF_DATA_DIR
global subBox
global MOUNTEDFLAG

systemPath='C:\Users\eladz\Documents\workProjects';
mountPath='\\wsl.localhost\Ubuntu-24.04\home\eladzing\sshfsMounts\vera';
veraPath='/vera/ptmp/gc/eladzing';


DRACO_SIMPATH = '/virgotng/universe/IllustrisTNG';
%/virgo/simulations/IllustrisTNG';%    '/ptmp/apillepi/IllustrisTNG';
%HOME_SIMPATH = '/home/zinger/sshfsMounts/draco/IllustrisTNG/simulationOutput';
HOME_SIMPATH = [ mountPath '\IllustrisTNG\simulationOutput'];

DRACO_PRINTOUT_DIR = [veraPath '/IllustrisTNG/printout'];
%/ptmp/eladzing/TNG/printout';
HOME_PRINTOUT_DIR = [ systemPath '/IllustrisTNG/printout'];

DRACO_MATFILE_DIR = [veraPath '/IllustrisTNG/matFiles'];
%'/u/eladzing/IllustrisTNG/matFiles';%/ptmp/eladzing/TNG/matFiles';
HOME_MATFILE_DIR = [ systemPath '/matlab_scripts/IllustrisTNG/matFiles'];

DRACO_FIG_DIR =[veraPath '/IllustrisTNG/printout/figFiles'];
%'/ptmp/eladzing/TNG/printout/figFiles';
HOME_FIG_DIR = [ systemPath '/IllustrisTNG/printout/figFiles'];

DRACO_myPOSTPROCESSING = [veraPath '/IllustrisTNG/ezPostProcessing'];
%'/ptmp/eladzing/TNG/ezpostprocessing';
HOME_myPOSTPROCESSING = [ mountPath '/IllustrisTNG/ezPostProcessing'];

DRACO_METHODS = [veraPath '/IllustrisTNG/TNG-Variations'];
HOME_METHODS = [mountPath '/IllustrisTNG/TNG-Variations'];

HOME_JF_DATA_DIR = [ systemPath '/IllustrisTNG/jellyFish/zooniverse_data'];
DRACO_JF_DATA_DIR = 'NOT DEFINED!';

DRACOFLAG=false;
doMount=false;

subbox='';

zred0=99;

%% parse arguments
i=1;
while(i<=length(varargin))
    if isstring(varargin{i})
        varargin{i}=char(varargin{i});
    end
    
    switch(lower(varargin{i}))
        case{'100','tng100',100,101,'100-1'}
            simName='L75n1820TNG';
            simDisplayName='TNG100';
            LBox=75e3; %in kpc/h
        case{'100lo','100low','100-2',102}
            simName='L75n910TNG';
            simDisplayName='TNG100-2';
            LBox=75e3; %in kpc/h
        case{'100lolo','100lower','100-3',103}
            simName='L75n455TNG';
            simDisplayName='TNG100-3';
            LBox=75e3; %in kpc/h
            
        case{'300','tng300',300}
            simName='L205n2500TNG';
            simDisplayName='TNG300';
            LBox=205e3; %in kpc/h
        case{'300lo','300low','300-2',302}
            simName='L205n1250TNG';
            simDisplayName='TNG300-2';
            LBox=205e3; %in kpc/h
        case{'300lolo','300lower','300-3',303}
            simName='L205n625TNG';
            simDisplayName='TNG300-3';
            LBox=205e3; %in kpc/h
            
        case{'50','tng50',50}
            simName='L35n2160TNG';
            simDisplayName='TNG50';
            LBox=35e3; %in kpc/h
        case{'50lo','50low','50-2',52}
            simName='L35n1080TNG';
            simDisplayName='TNG50-2';
            LBox=35e3; %in kpc/h
        case{'50lolo','50lower','50-3',53}
            simName='L35n540TNG';
            simDisplayName='TNG50-3';
            LBox=35e3; %in kpc/h
        case{'50lololo','50lowest','50-4',54}
            simName='L35n270TNG';
            simDisplayName='TNG50-4';
            LBox=35e3; %in kpc/h
         
        case{'clusters','tngclusters'}
            simName='L680n8192TNG';
            simDisplayName='TNGclust';
             LBox=680e3; %in kpc/h
                        
            DRACO_SIMPATH = '/virgotng/mpia/TNG-Cluster';%   
            HOME_SIMPATH = [ mountPath '\IllustrisTNG\TNG-Cluster'];
        
%         case{'clusters_dm','tngclusters_dm','clustersdm','tngclustersdm'}
%             simName='L680n2048TNG_DM';
%             simDisplayName='TNGclustDM';
%              LBox=680e3; %in kpc/h
%             
%             
%             DRACO_SIMPATH = '/virgotng/mpia/TNG-Cluster';%    '/ptmp/apillepi/IllustrisTNG';          
%             HOME_SIMPATH = [ mountPath '\IllustrisTNG\TNG-Cluster'];

        case{'illustris','illustris1','illustris-1','ill','ill1'}
            DRACO_SIMPATH = '/virgo/simulations/Illustris';%    '/ptmp/apillepi/IllustrisTNG';
            HOME_SIMPATH = '/home/zinger/sshfsMounts/draco/IllustrisTNG/Illustris';
            
            simName='Illustris-1';
            simDisplayName='Illustris-1';
            LBox=75e3; %in kpc/h
            
        case{'illustris2','illustris-2','ill2'}
            DRACO_SIMPATH = '/virgo/simulations/Illustris';%    '/ptmp/apillepi/IllustrisTNG';
            HOME_SIMPATH = '/home/zinger/sshfsMounts/draco/IllustrisTNG/Illustris';
            
            simName='Illustris-2';
            simDisplayName='Illustris-2';
            LBox=75e3; %in kpc/h
            
        case{'illustris3','illustris-3','ill3'}
            DRACO_SIMPATH = '/virgo/simulations/Illustris';%    '/ptmp/apillepi/IllustrisTNG';
            HOME_SIMPATH = '/home/zinger/sshfsMounts/draco/IllustrisTNG/Illustris';
            
            simName='Illustris-3';
            simDisplayName='Illustris-3';
            LBox=75e3; %in kpc/h
            
        case{'sub0','subbox0'}
            subbox='subbox0';
        case{'sub1','subbox1'}
            subbox='subbox1';
        case{'sub2','subbox2'}
            subbox='subbox2';
        case{'sub','subbox'}
            i=i+1;
            sb=varargin{i};
            if ~ischar(sb)
                sb=num2str(sb);
            end
            subbox=['subbox' sb];
        case{'method','methodbox','35','tng35',35}
            i=i+1;
            meth=varargin{i};
            if ~ischar(meth)
                error('SET_ENV: method code must be a string');
            end
            simName=['L25n512_' meth];
            simDisplayName=['TNG35_' meth];
            LBox=25e3; %in kpc/h
            
            DRACO_SIMPATH = DRACO_METHODS; %'/virgotng/universe/TNG-Variations'; %
            HOME_SIMPATH =  HOME_METHODS; %'/home/zinger/sshfsMounts/draco/IllustrisTNG/simulationMethods';
        case{'methodres'}
            i=i+1;
            res=varargin{i};
            if ~ischar(res)
                res=num2str(res);
            end
            k1=strfind(simName,'n');
            k2=strfind(simName,'_');
            
            simName(k1+1:k2-1)=res;
            
            k=strfind(simDisplayName,'_');
            
            
            simDisplayName=[ simDisplayName(1:k-1) '_n' res simDisplayName(k:end)];
            
            
        case 'draco'
            DRACOFLAG=true;
        case {'nomount','nomounting'}
            doMount=false;
        otherwise
            error('Unknown parameter in set_env: %s',varargin{i})
    end
    i=i+1;
end

%% Check to see if we are in draco

pp=pwd;
DRACOFLAG= strcmp(pp(1:6),'/vera/');



%% Set the environmmet for accssing IllustrisTNG data
if ~DRACOFLAG
    
    mPath='/home/eladzing/sshfsMounts/vera';
    
    %fprintf('*** Checking if Mounting simulation files *** \n');
    
    % check to see if mounted
    mountFlag=system(['wsl mountpoint -q ' mPath]); % 0 - is mounted, 1 - is not mounted
    MOUNTEDFLAG=~mountFlag;
    
    if MOUNTEDFLAG
        fprintf('*** remote filesystem on VERA mounted at: %s ***\n', mountPath);
    else
        warning('*** remote filesystem on VERA  - NOT MOUNTED!  ***\n');
    end
    
    %% this part is unneeded
    %     if doMount
    %
    %         % check to see if mounted
    %         mountFlag=system(['wsl mountpoint -q ' mPath]); % 0 - is mounted, 1 - is not mounted
    %         mountFlag=~mountFlag;
    %         maxIter=10;
    %         cnt=0;
    %         while ~mountFlag && cnt<=maxIter  % keep trying till you make it
    %             cnt=cnt+1;
    %             sysCommand=['wsl sshfs -p 2002 eladzing@localhost:. ' mPath ' -o follow_symlinks -o allow_other'];
    %             %'sshfs eladzing@draco.mpcdf.mpg.de:../../ptmp/apillepi/IllustrisTNG ' mPath ' -o follow_symlinks'];
    %             system(sysCommand);
    %
    %             mountFlag=system(['wsl mountpoint -q ' mPath]); % 0 - is mounted, 1 - is not mounted
    %             mountFlag=~mountFlag;
    %         end
    %         fprintf('*** Remote filesystem mounted on %s  *** \n',mountPath);
    %     else
    %         mountFlag=false;
end


%% set paths

if ~DRACOFLAG
    DEFAULT_PRINTOUT_DIR = HOME_PRINTOUT_DIR;
    DEFAULT_MOVIE_DIR = [ HOME_PRINTOUT_DIR '/movies'];
    DEFAULT_FIG_DIR = HOME_FIG_DIR;
    DEFAULT_MATFILE_DIR = HOME_MATFILE_DIR;
    DEFAULT_REMOTE_MATFILE_DIR = '/home/zinger/sshfsMounts/draco/IllustrisTNG/matFiles';
    myPOSTPROCESSING = HOME_myPOSTPROCESSING;
    JF_DATA_DIR = HOME_JF_DATA_DIR;
    SIMPATH = HOME_SIMPATH;
else
    DEFAULT_PRINTOUT_DIR = DRACO_PRINTOUT_DIR;
    DEFAULT_MOVIE_DIR = [ DRACO_PRINTOUT_DIR '/movies'];
    DEFAULT_FIG_DIR = DRACO_FIG_DIR;
    DEFAULT_MATFILE_DIR = DRACO_MATFILE_DIR;
    myPOSTPROCESSING = DRACO_myPOSTPROCESSING;
    JF_DATA_DIR = DRACO_JF_DATA_DIR;
    SIMPATH = DRACO_SIMPATH;
    mountFlag=true;
    
end

BASEPATH =[SIMPATH '/' simName '/output'];


POSTPROCESSING = [SIMPATH '/' simName '/postprocessing'];




% set cosmology
illustris.utils.set_cosmology(MOUNTEDFLAG);

%% snapshot - redshift info

load([DEFAULT_MATFILE_DIR '/snaps2redshift.mat']);

global fullSnapMask
global snapRedshifts

fullSnapMask=fullSnapMask0;
snapRedshifts=snapRedshifts0;

if contains(simDisplayName,'TNG35')
    zred0=4;
    snapRedshifts=[4 3 2 1 0];
    fullSnapMask=true(size(snapRedshifts));
end

% analyze sub-boxes
if ~isempty(subbox)
    BASEPATH=[BASEPATH '/' subbox];
    simDisplayName=[simDisplayName '-' subbox(1:3) subbox(end)];
    
    % enter the subbox sizes in units of ckpc/h
    switch simDisplayName
        case {'TNG100-sub0','TNG100-2-sub0','TNG100-3-sub0'}
            subBox.name='subbox0';
            subBox.center=1000.*[9 17 63];
            subBox.boxSize=7500;
            
        case {'TNG100-sub1','TNG100-2-sub1','TNG100-3-sub1'}
            subBox.name='subbox1';
            subBox.center=1000.*[37 43.5 67.5];
            subBox.boxSize=7500;
            subBox.Nsnaps=4380;
            
        case {'TNG300-sub0','TNG300-2-sub0','TNG300-3-sub0'}
            subBox.name='subbox0';
            subBox.center=1000.*[44 49 148];
            subBox.boxSize=15000;
            
        case {'TNG300-sub1','TNG300-2-sub1','TNG300-3-sub1'}
            subBox.name='subbox1';
            subBox.center=1000.*[20 175 15];
            subBox.boxSize=15000;
            
        case {'TNG300-sub2','TNG300-2-sub2','TNG300-3-sub2'}
            subBox.name='subbox1';
            subBox.center=1000.*[169 97.9 138];
            subBox.boxSize=10000;
            
            
        case {'TNG50-sub0','TNG50-2-sub0','TNG50-3-sub0','TNG50-4-sub0'}
            subBox.name='subbox0';
            subBox.center=1000.*[26 10 26.5];
            subBox.boxSize=4000;
            
        case {'TNG50-sub1','TNG50-2-sub1','TNG50-3-sub1','TNG50-4-sub1'}
            subBox.name='subbox1';
            subBox.center=1000.*[12.5 10 22.5] ;
            subBox.boxSize=4000;
            
        case {'TNG50-sub2','TNG50-2-sub2','TNG50-3-sub2','TNG50-4-sub2'}
            subBox.name='subbox1';
            subBox.center=1000.*[7.3 24.5 21.5];
            subBox.boxSize=5000;
            
        otherwise
            error('SET_ENV - sub box %s not implemented or non-existant',simDisplayName);
    end
    
    switch simDisplayName
        case {'TNG100-sub0','TNG100-sub1'}
            subBox.Nsnaps=7908;
            snapRedshifts=snapRedshifts_TNG100_sub;
        case {'TNG100-2-sub0','TNG100-2-sub1'}
            subBox.Nsnaps=4380;
        case {'TNG100-3-sub0','TNG100-3-sub1'}
            subBox.Nsnaps=2431;
            
        case {'TNG300-sub0','TNG300-sub1','TNG300-sub2'}
            subBox.Nsnaps=2449;
        case {'TNG300-2-sub0','TNG300-2-sub1','TNG300-2-sub2'}
            subBox.Nsnaps=3045;
        case {'TNG300-3-sub0','TNG300-3-sub1','TNG300-3-sub2'}
            subBox.Nsnaps=2050;
            
            %         case {'TNG50-sub0','TNG50-sub1','TNG50-sub2'}
            %             subBox.Nsnaps=;
            %         case {'TNG50-2-sub0','TNG50-2-sub1','TNG50-2-sub2'}
            %             subBox.Nsnaps=;
            %         case {'TNG50-3-sub0','TNG50-3-sub1','TNG50-3-sub2'}
            %             subBox.Nsnaps=;
            %         case {'TNG50-4-sub0','TNG50-4-sub1','TNG50-4-sub2'}
            %             subBox.Nsnaps=;
        otherwise
            error('SET_ENV - sub box %s not implemented or non-existant',simDisplayName);
            
    end
    
end

basePath=BASEPATH;



%% set some useful unit conversion factors
illustris.utils.set_illUnits(zred0) % set units to z=0
%
% global illUnits
% units;
% illUnits.massUnit=1e10/cosmoStruct.hub;  % Change simulation units to solar mass
% illUnits.lengthUnit=1/cosmoStruct.hub; % change units to kpc
% illUnits.volumeUnit=(illUnits.lengthUnit)^3; % change volume to kpc^3
% illUnits.densityUnit=illUnits.massUnit/illUnits.volumeUnit; % Change simulation units to solar mass / kpc^3
%
%
% % convert simulation density to number density in cm^-3
% illUnits.numberDensityFactor = (illUnits.densityUnit*Units.Ms/Units.kpc^3)/(cosmoStruct.muMass.*Units.mp);
%
% % convert density to Hydrogen number density in cm^-3
% illUnits.HydrogenNumberDensityFactor= (illUnits.densityUnit*Units.Ms/Units.kpc^3)*(cosmoStruct.Hfraction/Units.mp);
%
% % convert solarmass/kpc^3 to number density
% %illUnits.numberDensityFactor




fprintf('*** You are now analyzing the %s Simulation. ENJOY!  *** \n',simDisplayName);
end