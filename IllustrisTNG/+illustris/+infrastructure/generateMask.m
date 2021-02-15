function res = generateMask(varargin)
%GENERATEMASK generate a mask for 'good' galaxies

% defaults
global BASEPATH
global simDisplayName


centralFlag=false;
satFlag=false;
snap=99;
massThresh=10^9;
% haloMassThresh=10^11;

gasFlag=false;
dmFlag=false;

% parse argumnts
i=1;
while(i<=length(varargin))
    switch(lower(varargin{i}))
        case{'subs'}
            i=i+1;
            subs=varargin{i};
        case{'fofs'}
            i=i+1;
            fofs=varargin{i};
        case{'centrals'}
            centralFlag=true;
        case{'sats'}
            satFlag=true;
        case{'snap'}
            i=i+1;
            snap=varargin{i};
        case{'gas','hasgas'}
            gasFlag=true;
                    case{'dm','hasdm'}
            dmFlag=true;
        case{'mass','masshresh','thresh','threshold'}
            i=i+1;
            massThresh=varargin{i};
%         case{'halo','halomass'}
%             i=i+1;
%             haloMassThresh=varargin{i};
        otherwise
            error('generateMask - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end

%% set unit conversion
simUnits=illustris.utils.calc_illUnits(snap);

%% load subs & fofs
if ~exist('fofs','var')
    fofs=illustris.groupcat.loadHalos(BASEPATH,snap);
end

if ~exist('subs','var')
    subs=illustris.groupcat.loadSubhalos(BASEPATH,snap);
end
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);

%% define global mask

if ~contains(simDisplayName,'TNG35')
    
    baseMask=subs.SubhaloFlag & ...  % legal Subhalo
        subsInfo.DM10perc & ...      % avoid spurious clumps
        subsInfo.hasStars & ...      % galaxy which contians stars
        subsInfo.hostHasVirial;      % host has virial properties defined
else
    baseMask=subsInfo.DM10perc & ...      % avoid spurious clumps
        subsInfo.hasStars & ...      % galaxy which contians stars
        subsInfo.hostHasVirial;      % host has virial properties defined
end

%% define stellar mass mask
% if ~exist('massThresh','var')
%     switch simDisplayName
%         case'TNG100'
%             massThresh=10^9.5;
%         case 'TNG300'
%             massThresh=10^10.4;
%     end
% end

galMass= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*simUnits.massUnit; % stellar mass within 2*rhalf
massMask=galMass>massThresh;

res=baseMask & massMask;

% hostMass=fofs.Group_M_Crit200.*simUnits.massUnit;
% hostMask
%% gas mask 
if gasFlag
    res=res & subsInfo.hasGas;
end


if dmFlag
    res=res & subsInfo.hasDM;
end

%% centrals or satellite 
if centralFlag
    res=res & subsInfo.isCentral;
elseif satFlag
    res=res & ~subsInfo.isCentral; 
end



end

