function res=rscale_mmw_array(mstar,varargin)

%% align mstar correctly

sz=size(mstar);
if sz(1)==1
    mstar=mstar';
end

mask=true(size(mstar));
units;
maxIter=100;
GG=Units.G*Units.km^-2*Units.kpc^-1*Units.Ms;

% default values
fg=0.*ones(size(mstar)); %no gas in the disk
beta=1.*ones(size(mstar)); % gas scale radius same as stellar
xi=1000.*ones(size(mstar)); % bulge is a point mass
BT=0.*ones(size(mstar)); % no bulge component
%mdLimit=0.05; %maximal value of md
mdDefualt=0.03;
eps=1e-3;

%flags
lambdaFlag=false;
cvirFlag=false;
MvFlag=false;
mstarFlag=any(mstar>0);
mdFlag=false;
fbFlag=false;
btFlag=false;
jdFlag=false;
plotFlag=true;
converged=false;

i=1;
while i<=length(varargin)
    switch varargin{i}
        case 'fg'
            i=i+1;
            fg=varargin{i};
            fg=checkArgs(fg,mstar);
        case 'beta'
            i=i+1;
            beta=varargin{i};
            beta=checkArgs(beta,mstar);
        case {'BT','B/T'}
            i=i+1;
            BT=varargin{i};
            BT=checkArgs(BT,mstar);
            btFlag=true;
        case 'fb'
            i=i+1;
            fb=varargin{i};
            fb=checkArgs(fb,mstar);
            fbFlag=true;
        case 'xi'
            i=i+1;
            xi=varargin{i};
            xi=checkArgs(xi,mstar);
        case 'md'
            i=i+1;
            md=varargin{i};
            md=checkArgs(md,mstar);
            mdFlag=true;
            %         case {'md_limit','mdLimit'}
            %             i=i+1;
            %             mdLimit=varargin{i};
        case 'jd'
            i=i+1;
            jd0=varargin{i};
            jd0=checkArgs(jd0,mstar);
            jdFlag=true;
        case 'Mv'
            i=i+1;
            Mv=varargin{i};
            Mv=checkArgs(Mv,mstar);
            MvFlag=true;
        case {'cv','c','cvir'}
            i=i+1;
            c=varargin{i};
            c=checkArgs(c,mstar);
            cvirFlag=true;
        case 'lambda'
            i=i+1;
            lambda=varargin{i};
            lambda=checkArgs(lambda,mstar);
            lambdaFlag=true;
        case 'eps'
            i=i+1;
            eps=varargin{i};
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
    lambda=lambda_prime(mstar); % spin parameter
end


% find md, mstar Mv constellation
if mstarFlag
    if MvFlag
        md=mstar.*(1+fg)./Mv;
    elseif mdFlag
        Mv=mstar.*(1+fg)./md;
    else
        md=mdDefualt.*ones(size(mstar));
        %md=min(lambda,mdLimit);
        %baryon mass in disk (mstar+mgas)/Mv
        % assumed to be equal to spin parameter
        Mv=mstar.*(1+fg)./md;
    end
elseif MvFlag %mstar is not given by Mv is
    if mdFlag
        mstar=md.*Mv./(1+fg);
    else
        md=mdDefualt.*ones(size(mstar)); % min(lambda,mdLimit);
        mstar=md.*Mv./(1+fg);
    end
else
    error('rscale_mmw: neither Mv nor mstar has been given')
end

%% Halo info
if ~cvirFlag
    c=cvir_Mvir(Mv,0,'random'); % randomly select with scatter of 5% around formula (normal distribution)
end
[rvir,~,~,vvir]=calculate_virials('mv',Mv);
rvir=rvir.*1e3; % change to kpc
%A=Mv./(log(1+c)-c./(1+c));
%B=c./rvir;

% bulge stuff
if fbFlag && btFlag
    error('rscale_mmw: both fb and BT are given')
elseif fbFlag
    BT=fb./(fb+1);
end

Mb=BT.*mstar.*(1+fg)./(1-BT);
fb=Mb./mstar;
mb=Mb./Mv;

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
%lambdaLimit=max(0.1.*md,0.005.*ones(size(md)));
lambdaLimit=0.011; %max(0.1.*md,0.01.*ones(size(md)));
lambdaMask=lambda>lambdaLimit;
%lambdaLimit=0.005.*ones(size(md));
lambda=max(lambda,lambdaLimit);

fac=1;
rd_0=fac.*0.5.*sqrt(2).*lambda.*jd./md.*(1+fg).*rvir;


r0=logspace(-3.3,0.18,350); % 0.01:0.01:3;
rr=bsxfun(@times,r0,rvir);

rds=zeros(length(rd_0),maxIter); %zeros(length(numIter)+1,1);
rd=rd_0;
rds(:,1)=rd;
i=1;

%repmat
%AR=repmat(A,1,length(r0));
%BR=repmat(B,1,length(r0));
mdR=repmat(md,1,length(r0));
MbR=repmat(Mb,1,length(r0));
xiR=repmat(xi,1,length(r0));
MvR=repmat(Mv,1,length(r0));
cR=repmat(c,1,length(r0));
rvR=repmat(rvir,1,length(r0));
betaR=repmat(beta,1,length(r0));
fgR=repmat(fg,1,length(r0));
MsR=repmat(mstar,1,length(r0));

A=MvR./(log(1+cR)-cR./(1+cR));
B=cR./rvR;
mbR=MbR./MvR;

rdR=zeros(size(MsR));
ri=rdR;
r1=rdR;
r2=rdR;

vcdmSq=rdR;
vcbSq=rdR;
vcdkSq=rdR;
vc=rdR;
mDM=rdR;

b1=rdR;
bBeta=rdR;

eta=rdR;
integrand=rdR;

while i<maxIter
    i=i+1;
    maskR=repmat(mask,1,length(r0));
    rdR(mask,:)=repmat(rd(mask),1,length(r0));
    
%     if any(isnan(rdR))
%         pause
%     end
    
    %sigma=mstar./(2*pi*rd.^2);
    
    % solve adiabatic contraction
    ri(maskR)=arrayfun(@solve_MMW_NFW_one,rr(maskR),MvR(maskR),rvR(maskR),cR(maskR),MsR(maskR),rdR(maskR),fgR(maskR),betaR(maskR),MbR(maskR),xiR(maskR));
    %ri=solve_MMW_NFW_arr(rr,Mv,rvir,c,mstar,rd,fg,beta,Mb,xi);
    
    
    % find circular velocity
    
    vcdmSq(maskR)=GG.*((1-mdR(maskR)-mbR(maskR)).*A(maskR).*(log(1+B(maskR).*ri(maskR))-B(maskR).*ri(maskR)./(1+B(maskR).*ri(maskR))))./rr(maskR);
    mDM(maskR)=(1-mdR(maskR)-mbR(maskR)).*A(maskR).*(log(1+B(maskR).*ri(maskR))-B(maskR).*ri(maskR)./(1+B(maskR).*ri(maskR)));
    
    vcbSq(maskR)=(GG*MbR(maskR)).*rr(maskR)./(rdR(maskR)./xiR(maskR)+rr(maskR)).^2;
    
    r1(maskR)=0.5.*rr(maskR)./rdR(maskR);
    r2(maskR)=0.5.*rr(maskR)./rdR(maskR).*betaR(maskR);
    b1(maskR)=besseli(0,r1(maskR)).*besselk(0,r1(maskR))-besseli(1,r1(maskR)).*besselk(1,r1(maskR));
    bBeta(maskR)=besseli(0,r2(maskR)).*besselk(0,r2(maskR))-besseli(1,r2(maskR)).*besselk(1,r2(maskR));
    vcdkSq(maskR)=2.*GG.*MsR(maskR)./rdR(maskR).*r1(maskR).^2.*(b1(maskR)+fgR(maskR).*betaR(maskR).^3.*bBeta(maskR));
    %vcdkSq=GG.*pi.*sigma.*rd.*(rr./rd).^2.*(b1+fg.*beta.^3.*bBeta);
    
    vc(maskR)=sqrt(vcdmSq(maskR)+vcdkSq(maskR)+vcbSq(maskR));
    if any(isnan(vc))
        pause
    end
    % calculate integral
    eta(maskR)=2*r1(maskR);
    integrand(maskR)=eta(maskR).^2.*vc(maskR).*(exp(-1.*eta(maskR))+fgR(maskR).*betaR(maskR).^2.*exp(-1.*eta(maskR).*betaR(maskR)));
    integrand(integrand<1e-10)=0;
    
    zeta=zeros(size(mstar));
    
    for j=1:length(mstar)
        if mask(j)
            zeta(j)=trapz(eta(j,:),integrand(j,:))./vvir(j);
        end
    end
    
    
    rd(mask)=sqrt(2).*lambda(mask).*jd(mask)./md(mask).*(1+fg(mask)).*rvir(mask)./zeta(mask);
    rds(:,i)=rd;
    
%     if any(isnan(rd))
%         pause
%     end
    
    
    % check to leave loop
    converged=abs(diff(rds(:,i-1:i),1,2))./rd < eps;
    mask=~converged;
    
    if all(converged)
        break
    end
    
end

if any(~converged)
    fprintf(1,'Some values not converged');
    figure
    plot(rds(~converged,1:i)','.--');
    xlabelmine('iterations')
    ylabelmine('$R_d\,[\mathrm{kpc}]$')
    set(gca,'box','on','FontSize',12)
end

res.rd=rds(:,i);
res.Mv=Mv;
res.md=md;
res.BT=BT;
res.cv=c;
res.lambda=lambda;
res.converged=converged;
res.lambdaMask=lambdaMask;
res.mDM=mDM;
res.vc=vc;
res.rr=rr;
%res.eps=max(sqrt(vcdkSq))/sqrt(GG*mstar*(1+fg+fb)/rds(end));


if plotFlag
    figure
    plot(rds(:,1:i)','o--','linewidth',2);
    xlabelmine('iterations')
    ylabelmine('$R_s\,[\mathrm{kpc}]$')
    set(gca,'box','on','FontSize',14)
    grid
end
end

function res=checkArgs(arr,mstar)

if length(arr)==1
    arr=arr.*ones(size(mstar));
elseif size(arr')==size(mstar)
    arr=arr';
elseif size(arr)~=size(mstar)
    error('rscale_mmw_array: %s must be the same size as mstar',inputname(1));
end
res=arr;
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