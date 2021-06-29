%% batch Maps


doFlag=true;

while doFlag
    
    if ~exist('skipFlag','var')
        skipFlag=false;
    end
    if ~skipFlag
        %load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat');
        load('/home/zinger/workProjects/matlab_scripts/cluster/mat_files/shockedge.mat');
        list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
        
        %choose=chooseClusterProjection;
        cl=chooseCluster
        
        bmap='*RdYlBu';
        
        
        
        shk=shockedge_a1;
        shkM1=shockedgeMask1_a1;
        shkM2=shockedgeMask2_a1;
        
        
        boxx=8;
        new_env(cl);
        
        global NCELL
        global hub
        global zred
        global CLUSTER
        
        
        csSt.cluster=CLUSTER;
        k=strfind(CLUSTER,'CL')+2;
        clNum=str2num(CLUSTER(k:end));
        ind=find(list==clNum);
        
        c1=shk{ind,4};
        c2=shk{ind,5};
        circ=cat(2,c1(shkM1(ind,:)),c2(shkM2(ind,:)));
        
        
        %% perliminaries
        
        units;
        [xPos,yPos,zPos]=mk_position_cube(boxx,'hub',1.0);
        fluxnorm=(0.1117.*(get_mvir()./1e15).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
        gm=5/3;
        [Mg200 , ~, ~]=read_Mass_Profiles(get_rvir);
        [vx,vy,vz] = get_velocities(boxx);
        rcube2=xPos.^2+yPos.^2+zPos.^2 ; % r^2 cube in Mpc
        vr=(vx.*xPos+vy.*yPos+vz.*zPos)./sqrt(rcube2) ; %radial velocity
        
        dmro=smooth3(RHODM(boxx),'gaussian',7,3);
        ro=RHOG(boxx);
        tmp=T(boxx);
        ent=S(boxx).*Units.factors.f_ent;
        vtot=sqrt(vx.^2+vy.^2+vz.^2);
        %flux=ro.*vr.*rcube2./Mg200.*(km/Mpc*Gyr)./fluxnorm ;
        
        % boxx 8
        
        prj=chooseProjection;
    end
    
   
    
    csSt.prj=prj;
    
    
    
    %% plot reference plots
    thick=0.1;
    mkmap(boxx,'type','flux',prj,'clims',[-3 3],'vfield','dilute',6,'thick',thick,...
        'marks',circ(2),'brewer',bmap,'labels','full','title','no');
    hf(1)= gcf ;
    mkmap(boxx,'type','numberdensity',prj,'clims',[-9 -3],'vfield','dilute',6,'thick',thick,...
        'marks',circ(2),'brewer',bmap,'labels','full','title','no');
    hf(2)= gcf ;
    mkmap(boxx,'type','temp',prj,'clims',[4 8],'vfield','dilute',6,'thick',thick,...
        'marks',circ(2),'brewer',bmap,'labels','full','title','no');
    hf(3)= gcf ;
    mkmap(boxx,'type','ent',prj,'clims',[0 5],'vfield','dilute',6,'thick',thick,...
        'marks',circ(2),'brewer',bmap,'labels','full','title','no');
    hf(4)= gcf ;
    mkmap(boxx,'type','dm',prj,'clims',[7 15],'thick',thick,...
        'marks',circ(2),'brewer',bmap,'labels','full','title','no');
    hf(5)= gcf ;
    
    
    
    %% setup folders
    
    
    stat=true;
    pTag=0;
    while(stat)
        pTag=pTag+1;
        %printoutdir=sprintf('C:\\Users\\eladzing\\Documents\\cluster\\printout\\filaments\\crossSecs\\%s\\cs%s',CLUSTER,num2str(pTag));
        
        printoutdir=sprintf('/home/zinger/workProjects/cluster/printout/filaments/crossSecs/%s/cs%s',CLUSTER,num2str(pTag));
        
        [s,mess,messid]=mkdir(printoutdir);
        stat=strcmp(messid,'MATLAB:MKDIR:DirectoryExists');
    end
    csSt.printoutDir=printoutdir;
    csSt.pTag=pTag;
    %% draw line and get poition
    
    figure(hf(1));
    
    hl=imline;
    
    linePos=hl.getPosition;
    
    
    
    %% find cells closer to line than d.
    switch(lower(prj))
        case {'xy','yx'}
            p1=[linePos(1,1) linePos(1,2) 0];
            p2=[linePos(2,1) linePos(2,2) 0];
            ip=[1 2 3];
        case {'zy','yz'}
            p1=[0 linePos(1,2) linePos(1,1)];
            p2=[0 linePos(2,2) linePos(2,1)];
            ip=[3 2 1];
        case {'xz','zx'}
            p1=[linePos(1,1) 0 linePos(1,2)];
            p2=[linePos(2,1) 0 linePos(2,2)];
            ip=[1 3 2];
    end
    
    csSt.point1=p1;
    csSt.point2=p2;
    
    
    cellsize=boxx./hub./NCELL;
    dDist=1.5.*cellsize;
    
    csSt.dDist=dDist;
    
    dP=p2-p1;
    dPLength=sqrt(sum(dP.^2));
    xx1=(xPos-p1(1));xx2=(xPos-p2(1));
    yy1=(yPos-p1(2));yy2=(yPos-p2(2));
    zz1=(zPos-p1(3));zz2=(zPos-p2(3));
    
    csSt.dP=dP;
    
    tP=(xx1.*dP(1)+yy1.*dP(2)+zz1.*dP(3))./dPLength.^2;
    dd=sqrt((yy1.*zz2 - zz1.*yy2).^2 + ...
        (zz1.*xx2 - xx1.*zz2).^2 +...
        (xx1.*yy2 - yy1.*xx2).^2)./dPLength;
    
    mask=(tP>=0 & tP<=1) & dd<dDist;
    
    tt=tP(mask);
    csSt.tParam=tt;
    
    clear dd tP xx1 yy1 zz1 xx2 yy2 zz2
    
    
    hp=0.5.*(p1+p2);
    rh=hp./sqrt(sum(hp.^2));
    
    v1=dP;
    v2=v1;
    v2(ip(3))=v2(ip(3))+1;
    
    vn=cross(v1,v2);
    vn=vn./sqrt(sum(vn.^2));
    
    vn=-1.*sign(dot(vn,rh)).*vn;
    
    csSt.vNorm=vn;
    %plot mask
    
    pp1=p1+vn.*dDist;
    pp2=p2+vn.*dDist;
    pp3=p2-vn.*dDist;
    pp4=p1-vn.*dDist;
    
    name={'flux','rho','tmp','ent','dm'};
    for i=1:length(hf)
        figure(hf(i))
        patch([pp1(ip(1)) pp2(ip(1)) pp3(ip(1)) pp4(ip(1))],...
            [pp1(ip(2)) pp2(ip(2)) pp3(ip(2)) pp4(ip(2))],[1 0 0],...
            'FaceColor','none','edgeColor',[1 0 0],'linewidth',1,'linestyle','-')
        hold on
        quiver(p1(ip(1)),p1(ip(2)),dP(ip(1)),dP(ip(2)),0,'color','b','linewidth',1.1,'maxheadsize',0.6)
        fname=sprintf('%s_filCS_map_%s',CLUSTER,name{i});
        printout_fig(gcf,fname,'dir',printoutdir,'png');
    end
    %% get the necessary data
    
    
    % density
    roT=ro(mask).*(Units.Ms/Units.mm/Units.Mpc^3);
    
    % temperature
    tmpT=tmp(mask);
    
    % entropy
    entT=ent(mask); % in units of KeV cm^2;
    
    % pressure
    preT=(gm*Units.kb/Units.mm).*roT.*(Units.Ms/Units.Mpc^3).*tmpT;
    
    % flux
    fluxT=roT.*vr(mask).*rcube2(mask);  %flux(mask);
    
    % total velocity
    vt=vtot(mask)./get_vvir;
    
    % local Mach
    cso=(gm*Units.kb/Units.mm.*tmpT).^0.5./Units.km;  % ./1e5;  % so
    mac=vtot(mask)./cso;
    
    % perpindicular velocity
    
    vp=(vx.*vn(1)+vy.*vn(2)+vz.*vn(3))./get_vvir;
    vp=vp(mask);
    
    vrT=vr(mask)./get_vvir;
    
    mu=(vx.*vn(1)+vy.*vn(2)+vz.*vn(3))./sqrt(vx.^2+vy.^2+vz.^2);
    muT=mu(mask);
    % dm
    dmT=dmro(mask);
    
    
    csSt.ro=roT;
    csSt.tmp=tmpT;
    csSt.ent=entT;
    csSt.pres=preT;
    csSt.flux=fluxT;
    csSt.mach=mac;
    csSt.vt=vt;
    csSt.vr=vrT;
    csSt.vp=vp;
    csSt.mu=muT;
    csSt.dm=dmT;
    
    dPl=ceil(1e3.*dPLength./hub);
    csSt.dPLen=dPl;
    %% plot
    
    dPl=ceil(1e3.*dPLength./hub);
    csSt.dPLen=dPl;
    
    plot_filament_crossSection( csSt,'all' )
    
    save(sprintf('%s/cs.mat',printoutdir),'csSt');
    
    
    %% ask about the future
    answer = questdlg('Continue ?', ...
        'What Next?','Yes','No','No');
    switch answer
        case 'Yes'
            doFlag=true;
        case 'No'
            doFlag=false;
    end
    
    if doFlag
        
        answer = questdlg('Try Another Cross-Section ?', ...
            'What Next?', ...
            'Same Projection','Same Cluster, Different Projection','Different Cluster','Same Projection');
        
        switch answer
            case 'Same Projection'
                skipFlag=true;
                
            case 'Same Cluster, Different Projection'
                skipFlag=true;
                prj=chooseProjection;
            case 'Different Cluster'
                skipFlag=false;
            otherwise
                doFlag=false;
        end
        
        close all
        
    end
    
end




%
%
% xLab=sprintf('$t_p\\,[%s\\,\\mathrm{kpc}]$',num2str(dPl));
% % density
% figure
% semilogy(tt,roT,'.')
% xlabelmine(xLab)
% ylabelmine('$\log n_{gas}\,[\mathrm{cm^{-3}}]$');
% grid
% fname=sprintf('%s_filCS_rho_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % temperaturem
% semilogy(tt,tmpT,'.')
% xlabelmine(xLab)
% ylabelmine('$\log T\,[\mathrm{K}]$');
% grid
% fname=sprintf('%s_filCS_tmp_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % entropy
% figure
% semilogy(tt,entT,'.')
% xlabelmine(xLab)
% ylabelmine('$\log S\,[\mathrm{KeV\,cm^2}]$');
% grid
% fname=sprintf('%s_filCS_ent_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % pressure
% figure
% semilogy(tt,preT,'.')
% xlabelmine(xLab)
% ylabelmine('$\log P\,[\mathrm{dym\,cm^{-2}}]$');
% grid
% fname=sprintf('%s_filCS_press_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % mach
% figure
% semilogy(tt,mac,'.')
% xlabelmine(xLab)
% ylabelmine('local $\mathcal{M}$');
% grid
% fname=sprintf('%s_filCS_mach_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % vtot
% figure
% plot(tt,vt,'.')
% xlabelmine(xLab)
% ylabelmine('$v/V_{\mathrm{vir}}$');
% grid
% fname=sprintf('%s_filCS_vtot_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % vr
% figure
% plot(tt,vrT,'.')
% xlabelmine(xLab)
% ylabelmine('$v_r/V_{\mathrm{vir}}$');
% grid
% fname=sprintf('%s_filCS_vr_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % mu
% figure
% plot(tt,muT,'.')
% xlabelmine(xLab)
% ylabelmine('$\mu$');
% grid
% fname=sprintf('%s_filCS_mu_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % flux
% figure
% plot(tt,fluxT,'.')
% xlabelmine(xLab)
% ylabelmine('$\dot{M}_{gas}/\mathrm{d}\Omega$ [arbit. units]');
% grid
% fname=sprintf('%s_filCS_flux_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % perpindicular velocity
% figure
% plot(tt,vp,'.')
% xlabelmine(xLab)
% ylabelmine('$v_{\perp}/v_{\mathrm{vir}}$');
% grid
% fname=sprintf('%s_filCS_vperp_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);
%
% % dm
% figure
% semilogy(tt,dmT,'.')
% xlabelmine(xLab)
% ylabelmine('$\rho_{\mathrm{DM}}\,[\mathrm{M_\odot \, Mpc^{-3}}]$')
% grid
% fname=sprintf('%s_filCS_dm_cs%s',CLUSTER,num2str(pTag));
% printout_fig(gcf,fname,'dir',printoutdir);

