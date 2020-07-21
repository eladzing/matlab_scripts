%% defining a Hernquist halo object.
% all radial positions are assumed to be in units of the virial radius

classdef HERNQUIST
    properties
        Mh
        Rs
        GG=4.2997e-06; % gravitational constant in Msun kpc km/sec
        header='Hernquist model in units of solarMass, Kpc, km/sec';
        
    end
    methods
        % constructor
        function obj=HERNQUIST(varargin)
        
        mv=[];
        %cc=[];
        %rv=[];
        rs=[];
        %cosmoSt=[];
        
        %delv=200;
        %rhoType='crit';
        %zred=0;
        
        i=1;
        while i<=length(varargin)
            switch(lower(varargin{i}))
                case{'mv','mvir','m','mh','mb'}
                    i=i+1;
                    mv=varargin{i};
%                 case{'cc','cvir','c','cv'}
%                     i=i+1;
%                     cc=varargin{i};
%                     
%                 case{'rv','rvir','r'}
%                     i=i+1;
%                     rv=varargin{i};
                  case{'rs','rscale','rb'}
                    i=i+1;
                    rs=varargin{i};
%                                        
%                 case {'delv','delta','delta_vir','deltavir'}
%                     i=i+1;
%                     delv=varargin{i};
%                     
%                     
%                 case {'mean','rho_mean','rhomean'}
%                     rhoType='mean';
%                 case {'crit','rho_crit','rhocrit'}
%                     rhoType='crit';
%                     
%                 case {'cosmostruct','cosmo'}
%                     i=i+1;
%                     cosmoSt=varargin{i};
%                     if ~isstruct(cosmoSt)
%                         error('HERNQUIST: cosmology structure must be a structure')
%                     end
%                     
                otherwise
                    error('HERNQUIST: illegal argument: %s',varargin{i})
            end
            i=i+1;
        end
        
        if isempty(mv) ||  isempty(rs)
            error('HERNQUIST: must enter both  mass and scale radius ');
        end
        
                
%         if isempty(cosmoSt)
%             cosmoSt.Omm=0.3;
%             cosmoSt.Oml=0.7;
%             cosmoSt.hub=0.7;
%             cosmoSt.muMass=0.5926;
%         end
%         
        
        
%         if isempty(rv)
%             [rv, ~, tv, vv]=calculate_virials('mvir',mv,...
%                 'cosmoStruct',cosmoSt,'delv',delv','zred',zred,rhoType);
%             
%         else
%             [~, ~, tv, vv]=calculate_virials('mvir',mv,...
%                 'cosmoStruct',cosmoSt,'delv',delv','zred',zred,rhoType);
%         end
%         
%         if isempty(rs)
%             rs=rv/cc;
%         elseif isempty(cc)
%             cc=rv/rs;
%         end
%         
        
        obj.Mh=mv;  %in Msun
%         obj.Rvir=1000.*rv; %in kpc
%         obj.Tvir=tv; %in km/sec
%         obj.Vvir=vv; % in kelvin 
%         obj.cvir=cc; % unitless
        obj.Rs=rs;
        end
        
        
        % enclosed mass profile
        function res=mass(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rs;
                case 'mpc'
                    rr=rr./(obj.Rs/1000);
            end
        end
        
                
        res=obj.Mh.*rr.^2./(1+rr).^2;
                
        end
        
        % density profile
        function res=density(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rs;
                case 'mpc'
                    rr=rr./(obj.Rs/1000);
            end
        end
        
        
        fac=(obj.Mh/(2*pi*obj.Rs^3));
        res=fac./(rr.*(1+rr).^3);
        
                
        
        end
        
        function res=rho(obj,rr,runit)
        res=density(obj,rr,runit);
        end
        
                
        % potential
        function res=potential(obj,rr,runit)
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rs;
                case 'mpc'
                    rr=rr./(obj.Rs/1000);
            end
        end
        
                
        res=(-1*obj.GG*obj.Mh/obj.Rs)./(1+rr);
        
        end
        
        % circular velocity 
        function res=vcirc(obj,rr,runit)
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rs;
                case 'mpc'
                    rr=rr./(obj.Rs/1000);
            end
        end
        
        
        res=sqrt((obj.GG*obj.Mh/obj.Rs).*rr./(1+rr).^2);
        
        end
        
        
        
        % gravitaional accelaration
        function res=ggrav(obj,rr,runit)
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rs;
                case 'mpc'
                    rr=rr./(obj.Rs/1000);
            end
        end
        
        
        res=(-1*obj.GG*obj.Mh/obj.Rs^2)./(1+rr).^2;
        
        end
 
 
        
        
        
        
    end
end
