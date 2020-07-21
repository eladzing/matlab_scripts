%% defining a full galaxy model

classdef GALAXY
    properties
        StellarDisk
        GasDisk
        Bulge
        Halo
        Ms=[];
        Rd=[];
        fgs=[];
        beta=[];
        fbs=[];
        xi=[];
        Mh=[];
        cv=[];
        fg=[];
        gasMass
        rr=[];
    end
    
    methods
        %constructor
        function obj=GALAXY(varargin)
        Ms=[];
        Rd=[];
        fgs=[];
        beta=[];
        fbs=[];
        xi=[];
        Mh=[];
        cv=[];
        fg=[];
        rr=0:0.01:15;
        %parameter list: Ms, Rd, fgs, beta, fbs, xi, Mh, cv, fg
        %                disk    gas disk     bukge   halo
        
        % here comes a whole complicated thing about reading in all the parameters
        i=1;
        while i<=length(varargin)
            switch(lower(varargin{i}))
                case{'ms'}
                    i=i+1;
                    Ms=varargin{i};
                case{'rd'}
                    i=i+1;
                    Rd=varargin{i};
                case{'fgs'}
                    i=i+1;
                    fgs=varargin{i};
                case{'beta'}
                    i=i+1;
                    beta=varargin{i};
                case{'fbs'}
                    i=i+1;
                    fbs=varargin{i};
                case{'xi'}
                    i=i+1;
                    xi=varargin{i};
                case{'mh'}
                    i=i+1;
                    Mh=varargin{i};
                case{'cv'}
                    i=i+1;
                    cv=varargin{i};
                    
                case{'fg'}
                    i=i+1;
                    fg=varargin{i};
                case{'rr'}
                    i=i+1;
                    rr=varargin{i};
                otherwise
                    error('GALAXY - Illegal argument for this class: %s',varargin{i});
            end
            i=i+1;
        end
        
        % here we randomly select those parameters that were not given
        
        % build the different components
        obj.StellarDisk=EXPDISK('Ms',Ms,'Rd',Rd);
        
        Mg=fgs.*Ms;
        Rg=Rd./beta;
        obj.GasDisk=EXPDISK('Ms',Mg,'Rd',Rg);
        
        Mb=fbs.*Ms;
        Rb=Rd./xi;
        obj.Bulge=HERNQUIST('Mh',Mb,'Rs',Rb);
        
        if ~isempty(Mh)
            obj.Halo=NFW('Mv',Mh,'cv',cv);
        end
        
        % create gas mass profile
        obj.gasMass.rr=rr.*Rd; %in kpc
        obj.gasMass.mass=diff(mass(obj.GasDisk,obj.gasMass.rr,'kpc'));
        %obj.gasMass.sigma=GALAXY.calcSigma(obj.gasMass.rr,obj.gasMass.mass);
        
        obj.Ms=Ms;
        obj.Rd=Rd;
        obj.fgs=fgs;
        obj.beta=beta;
        obj.fbs=fbs;
        obj.xi=xi;
        obj.Mh=Mh;
        obj.cv=cv;
        obj.fg=fg;
       
        end
        
        
        %% gravitational acceleration
        function res=ggrav(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rd;
                case 'mpc'
                    rr=rr./(obj.Rd/1000);
            end
        end
        
        % radius is now in units of the stellar scale radius.
        % we recast it in kpc and let each component calculate
        
        rr=rr.*obj.Rd;
        
        res=ggrav(obj.StellarDisk,rr,'kpc')+ggrav(obj.GasDisk,rr,'kpc')+...
            ggrav(obj.Bulge,rr,'kpc')+ggrav(obj.Halo,rr,'kpc');
        
        end
        
        
        %% calculate binding force acting on gas disk
        function res=gasBindingForce(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rd;
                case 'mpc'
                    rr=rr./(obj.Rd/1000);
            end
        end
        
        rr=rr.*obj.Rd; % now in kpc
        
        res=ggrav(obj,rr,'kpc').*density(obj.GasDisk,rr,'kpc');
        
        
        end
        
        %% circular velocity
        function res=vcirc(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rd;
                case 'mpc'
                    rr=rr./(obj.Rd/1000);
            end
        end
        
        rr=rr.*obj.Rd;  % now in kpc
        
        res=vcirc(obj.StellarDisk,rr,'kpc')+vcirc(obj.GasDisk,rr,'kpc')+...
            vcirc(obj.Bulge,rr,'kpc')+vcirc(obj.Halo,rr,'kpc');
        
        
        
        end
        
        %% total enclosed mass
        %function res=enclosedMass(obj,rr)
        %end
    end
    methods (Static)
        
        %% calculate the total mass of stars formed per unit time baseb on kennicutt-schmidt
        function res=starFormationRate(rr,mass,varargin)
        
        randFlag=false;
        Afac=2.5;
        alfa=1.4;
        
        AfacScat=0.7;
        alfaScat=0.15;
        
        i=1;
        while i<=length(varargin)
            switch lower(varargin{i})
                case{'alf','alfa','alpha','alph'}
                    i=i+1;
                    alfa=varargin{i};
                case {'a','afac','fac'}
                    i=i+1;
                    Afac=varargin{i};
                case{'alfscatter','alfascatter','alphascatter','alphscatter'}
                    i=i+1;
                    alfaScat=varargin{i};
                case {'ascat','ascatter','afacscat','afacscatter','facscat','facscatter'}
                    i=i+1;
                    AfacScat=varargin{i};
                case {'rand','random','scatter'}
                    randFlag=true;
                otherwise
                error('STARFORMATIONRATE - Illegal argument %s',varargin{i});
            end
            i=i+1;
        end
        
        if randFlag
            Afac=Afac+AfacScat.*(2.*rand()-1);
            alfa=alfa+alfaScat.*(2.*rand()-1);
        end
            
        sigma=GALAXY.calcSigma(rr,mass);
        
        res=(Afac.*1e-4.*(sigma.*1e-6).^alfa).*GALAXY.dArea(rr); %in units of Msun/yr
        
        end
                
        %% calculate surface density from gas mass profile. 
        % mass profile is given as M(<r). 
        function res=calcSigma(rr,mm)
        
        res=mm./GALAXY.dArea(rr);
        
        end
        
        %% calculate ds of a ring
        function res=dArea(rr)
        res=pi.*(rr(2:end).^2-rr(1:end-1).^2);
        end
    end
end
