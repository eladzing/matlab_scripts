function res=rscale_mmw(mstar,varargin)

units;
numIter=100;
GG=G*km^-2*kpc^-1*Ms;

% default values
fg=0; %no gas in the disk
beta=1; % gas scale radius same as stellar
xi=1000; % bulge is a point mass
BT=0; % no bulge component
mdLimit=0.05; %maximal value of md
mdDefualt=0.1;

%flags
lambdaFlag=false;
cvirFlag=false;
MvFlag=false;
mstarFlag=mstar>0;
mdFlag=false;
fbFlag=false;
btflag=false;
jdFlag=false;
plotFlag=true;
converged=false;

i=1;
while i<=length(varargin)
    switch varargin{i}
        case 'fg'
            i=i+1;
            fg=varargin{i};
        case 'beta'
            i=i+1;
            beta=varargin{i};
        case {'BT','B/T'}
            i=i+1;
            BT=varargin{i};
            btFlag=true;
        case 'fb'
            i=i+1;
            fb=varargin{i};
            fbFlag=true;
        case 'xi'
            i=i+1;
            xi=varargin{i};
        case 'md'
            i=i+1;
            md=varargin{i};
            mdFlag=true;
        case {'md_limit','mdLimit'}
            i=i+1;
            mdLimit=varargin{i};
            
        case 'jd'
            i=i+1;
            jd0=varargin{i};
            jdFlag=true;
        case 'Mv'
            i=i+1;
            Mv=varargin{i};
            MvFlag=true;
        case {'cv','c','cvir'}
            i=i+1;
            c=varargin{i};
            cvirFlag=true;
        case 'lambda'
            i=i+1;
            lambda=varargin{i};
            lambdaFlag=true;
        case{'noplot','noshow'}
            plotFlag=false;
        case{'plot','show'}
            plotFlag=true;
        otherwise
            error('rscale_mmw: Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

%spin parameter
if ~lambdaFlag
    lambda=lambda_prime(1); % spin parameter
end


% find md, mstar Mv constellation
if mstarFlag
    if MvFlag
        md=mstar*(1+fg)/Mv;
    elseif mdFlag
        Mv=mstar*(1+fg)/md;
    else
        md=mdDefualt;
        %md=min(lambda,mdLimit);
        %baryon mass in disk (mstar+mgas)/Mv
        % assumed to be equal to spin parameter
        Mv=mstar*(1+fg)/md;
    end
elseif MvFlag %mstar is not given by Mv is
    if mdFlag
        mstar=md*Mv/(1+fg);
    else
        md=mdDefualt; % min(lambda,mdLimit);
        mstar=md*Mv/(1+fg);
    end
else
    error('rscale_mmw: neither Mv nor mstar has been given')
    
end

%% Halo info
if ~cvirFlag
    c=cvir_Mvir(Mv,0);
end
[rvir,~,~,vvir]=calculate_virials('mv',Mv);
rvir=rvir*1e3; % change to kpc
A=Mv/(log(1+c)-c/(1+c));
B=c/rvir;

% bulge stuff
if fbFlag && btFlag
    error('rscale_mmw: both fb and BT are given')
elseif fbFlag
    BT=fb/(fb+1);
end

Mb=BT*mstar*(1+fg)/(1-BT);
fb=Mb/mstar;
mb=Mb/Mv;

% jd models
if jdFlag
    if ischar(jd0)
        switch jd0
            case{'jdmd','md','jd=md','op1','def','default'}
                jd=md; %fraction of AM in disk
            case{'mdmb','md+mb','jd=md+mb','op2'}
                jd=md+mb;
            otherwise
                erorr('rscale_mmw: Illegal jd option')
        end
    elseif isnumeric(jd0)
        jd=jd0;
    else
        error('rscale_mmw: what did you put in jd?')
    end
else
    jd=md; %fraction of AM in disk
end


%% begin  iterative calculation
%needed for convergence
lambda=max(lambda,max(0.1*md,0.005));

fac=1;
rd_0=fac*0.5*sqrt(2)*lambda*jd/md*(1+fg)*rvir;





r0=0.01:0.01:3;
rr=r0.*rvir;
rds=[]; %zeros(length(numIter)+1,1);
rd=rd_0;
rds(1)=rd;
i=1;
while i<numIter
    
    %mdisk=mstar.*(1-exp(-1.*rr./rd).*(1+rr./rd)+fg.*(1-exp(-1.*rr./rd.*beta).*(1+rr./rd.*beta)));
    sigma=mstar/(2*pi*rd^2);
    
    % solve adiabatic contraction
  %  rrComp=rr(~isreal(rr));
  %  if ~isempty(rrComp)
  %      pause
  %  end

  
    ri=solve_MMW_NFW(rr,Mv,rvir,c,mstar,rd,fg,beta,Mb,xi);
    
    
    % find circular velocity
    vcdmSq=GG.*((1-md-mb).*A.*(log(1+B.*ri)-B.*ri./(1+B.*ri)))./rr;
    vcbSq=(GG*Mb).*rr./(rd/xi+rr).^2;
    
    b1=besseli(0,0.5.*rr./rd).*besselk(0,0.5.*rr./rd)-besseli(1,0.5.*rr./rd).*besselk(1,0.5.*rr./rd);
    bBeta=besseli(0,0.5.*rr./rd.*beta).*besselk(0,0.5.*rr./rd.*beta)...
        -besseli(1,0.5.*rr./rd.*beta).*besselk(1,0.5.*rr./rd.*beta);
    vcdkSq=GG.*pi.*sigma.*rd.*(rr./rd).^2.*(b1+fg.*beta^3.*bBeta);
    
    vc=sqrt(vcdmSq+vcdkSq+vcbSq);
    
    % calculate integral
    eta=rr./rd;
    integrand=eta.^2.*vc.*(exp(-1.*eta)+fg.*beta.^2.*exp(-1.*eta.*beta));
    
    zeta=trapz(eta,integrand)./vvir;
    
    rd=sqrt(2)*lambda*jd/md*(1+fg)*rvir./zeta;
    rds(end+1)=rd;
    
    % check to leave loop
    
    if abs(diff(rds(end-1:end)))/rd < 5e-4  %
        converged=true;
        break
    end
    i=i+1;
end

res.rd=rds(end);
res.Mv=Mv;
res.md=md;
res.BT=BT;
res.cv=c;
res.lambda=lambda;
res.converged=converged;
%res.eps=max(sqrt(vcdkSq))/sqrt(GG*mstar*(1+fg+fb)/rds(end));


if plotFlag
    figure
    plot(rds,'.--');
    xlabelmine('iterations')
    ylabelmine('$R_d\,[\mathrm{kpc}]$')
    set(gca,'box','on','FontSize',12)
    
end


% figure
% h=[];
% h(1)=plot(rr,sqrt(vcdmSq),'-r','DisplayName','DM','LineWidth',2);
% hold on
% h(2)=plot(rr,sqrt(vcdkSq),'-b','DisplayName','Disk','LineWidth',2);
% h(3)=plot(rr,sqrt(vcbSq),'-g','DisplayName','Bulge','LineWidth',2);
% h(4)=plot(rr,vc,'-k','DisplayName','Total','LineWidth',2);
% h(5)=plot(rd.*[1 1],ylim,':k','DisplayName','$R_d$','LineWidth',2.5);
% h(6)=plot(rvir.*[1 1],ylim,'--k','DisplayName','$R_{\mathrm{vir}}$','LineWidth',2);
% h(7)=plot(xlim,vvir.*[1 1],'--g','DisplayName','$V_{\mathrm{vir}}$','LineWidth',2);
% %ylim([10 200]);
% hl=legend(h);
% set(hl','Fontsize',14,'Interpreter','latex');
% xlabelmine('$r\,[\mathrm{kpc}]$');
% ylabelmine('$V_c\,[\mathrm{km/sec}]$');
% grid
% set(gca,'box','on','FontSize',12)
%
% figure
% h=[];
% h(1)=plot(rr,sqrt(vcdmSq),'-r','DisplayName','DM','LineWidth',2);
% hold on
% h(2)=plot(rr,sqrt(vcdkSq),'-b','DisplayName','Disk','LineWidth',2);
% h(3)=plot(rr,sqrt(vcbSq),'-g','DisplayName','Bulge','LineWidth',2);
% h(4)=plot(rr,vc,'-k','DisplayName','Total','LineWidth',2);
% h(5)=plot(rd.*[1 1],ylim,':k','DisplayName','$R_d$','LineWidth',2.5);
% h(6)=plot(rvir.*[1 1],ylim,'--k','DisplayName','$R_{\mathrm{vir}}$','LineWidth',2);
% h(7)=plot(xlim,vvir.*[1 1],'--g','DisplayName','$V_{\mathrm{vir}}$','LineWidth',2);
% xlim([0 25]);
% hl=legend(h);
% set(hl','Fontsize',14,'Interpreter','latex');
% xlabelmine('$r\,[\mathrm{kpc}]$');
% ylabelmine('$V_c\,[\mathrm{km/sec}]$');
% grid
% set(gca,'box','on','FontSize',12)