if perlimFlag
    global illUnits
    %global cosmoStruct
    units
      
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
    %load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/fofs_subs_TNG100_z0.mat')
    
    subsInfo=illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    
    mask=illustris.utils.generateMask('subs',subs,'fofs',fofs,'centrals');
    
    rr=(0.01:0.01:2);
end
tffProfile.polyfit.a=zeros(size(mask));
tffProfile.polyfit.b=zeros(size(mask));
tffProfile.polyfit.c=zeros(size(mask));
tffProfile.polyfit.Rsqr=zeros(2,length(mask));
tffProfile.polyfit.residuals=nan(6,length(mask));

% tffProfile.NFW.mv=zeros(size(mask));
% tffProfile.NFW.cv=zeros(size(mask));
% tffProfile.NFW.residuals=nan(6,length(mask));

fprintf('[');
stp=5;
cnt=10;
len=length(mask);

for i=1:len
    
    if ceil(100*(i./len))>=cnt
        fprintf('=');
        cnt=cnt+stp;
    end
    
     
    
    if ~mask(i)
        continue
    end
    
    idGal=i-1;
    idHost=subsInfo.hostFof(idGal+1);
    
    m=double([subs.SubhaloMassInRad(idGal+1) ...
        subs.SubhaloMassInMaxRad(idGal+1) ...
        subs.SubhaloMass(idGal+1)./2 ...
        fofs.Group_M_Crit500(idHost+1) ...
        fofs.Group_M_Crit200(idHost+1) ...
        fofs.Group_M_Mean200(idHost+1) ...
        ].*illUnits.massUnit);
    r=double([2*subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,idGal+1) ...
        subs.SubhaloVmaxRad(idGal+1) ...
        subs.SubhaloHalfmassRad(idGal+1) ...
        fofs.Group_R_Crit500(idHost+1) ...
        fofs.Group_R_Crit200(idHost+1) ...
        fofs.Group_R_Mean200(idHost+1) ...
        ].*illUnits.lengthUnit);
    
    wts=ones(size(m));
    wts(2)=0.5;
    wts(3)=0.2;
    
    mv=m(5);  % in Ms
    rv=r(5); % in kpc
    
    rhom=(3*m./(4*pi.*r.^3)).*Units.Ms./Units.kpc^3;
    tff0=sqrt(3*pi./(32*Units.G.*rhom))./Units.Gyr;
    
    
    %% fit 1 - polynomial
    % Set up fittype and options.
    
    xData=log10(r./rv)';
    yData=log10(tff0)';
    
    ft = fittype( 'poly2' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'Bisquare';
    opts.Weights = wts;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    tffProfile.polyfit.a(i)=fitresult.p1;
    tffProfile.polyfit.b(i)=fitresult.p2;
    tffProfile.polyfit.c(i)=fitresult.p3;
    
    tffProfile.polyfit.Rsqr(1,i)=gof.rsquare;
    tffProfile.polyfit.Rsqr(2,i)=gof.adjrsquare;
    
    
    tff1=10.^(fitresult.p1.*log10(r./rv).^2 + fitresult.p2.*log10(r./rv) + fitresult.p3);
    
    tffProfile.polyfit.residuals(:,i)=100.*(tff1./tff0-1)';
    
    
    %% fit 2 - NFW
    
    
%     cv=cvir_Mvir_200(mv,0,'hub',cosmoStruct.hub);
%     
%     halo=NFW('mv',mv,'cv',cv,'cosmo',cosmoStruct);
%     
%     %tff2=freefallTime(halo,rr);
%     
%     tff2=interp1(rr,freefallTime(halo,rr),r./rv);
%     tffProfile.NFW.mv(i)=mv;
%     tffProfile.NFW.cv(i)=cv;
%     tffProfile.NFW.residuals(:,i)=100.*(tff2./tff0-1)';
    
end

fprintf(']\n');


global DEFAULT_MATFILE_DIR
global simDisplayName
    
fname=sprintf('freeFallTime_profiles_snp%s_%s.mat',num2str(snap),simDisplayName);
save([DEFAULT_MATFILE_DIR '/' fname],'tffProfile','mask','-v7.3')
    

