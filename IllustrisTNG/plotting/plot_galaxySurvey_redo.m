%% plot images for galaxies
simName='100';
snap=99;
bp=illustris.set_env(simName);

global simDisplayName

if perlimFlag
    
    if readFlag
        fprintf(' *** Reading data *** \n');
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
        
        
    end
    
    
    massThresh=10^9; % threshold for *stellar* mass
    
    
    ssfr=illustris.utils.calc_ssfr(subs);
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    
    global illUnits
    
    massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
    
    
    
    %galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals','hasGas');
    
end


m9=[572671, 585431, 586813, 587600, 588650, 594081, 597343, 600544, 601580, 605760, 606203,...
    607550, 609036, 612501, 613305, 614090, 616295, 616552, 619292, 624453, 626573, 627692, ...
    632055, 633663, 635097, 639317, 640190, 642323, 645197, 647808, 649327, 653913, 654134, ...
    655967, 658143, 658517, 659027, 666746, 669775, 671560, 672233, 674513, 678411, 679057, ...
    684753, 687261, 702031];

m95=[465842, 540200, 546567, 551333, 552112, 554528, 554728, 560367, 565094, 570763, 572092,...
    579612, 585524, 587022, 592594, 593480, 595189, 596464, 597550, 599248, 599460, 605609,...
    606064, 606889, 608016, 610747, 617661, 619048, 622533, 632020, 637717];

m10=[444251, 455730, 459169, 460526, 470808, 471901, 488174, 488841, 490053, 497088, 499161,...
    507157, 511214, 511787, 516188, 517336, 518614, 522484, 522749, 523157, 530514, 532912,...
    534514, 535117, 540318, 540409, 541808, 542044, 556458, 566543];


m105=[325335, 350183, 372778, 387726, 390859, 399679, 401670, 405568, 414759, 421134, 426575,...
    431751, 434152, 438038, 441397, 444997, 447153, 452651, 452737, 456634, 462077, 467798,...
    469803, 470617, 472174, 473166, 479264, 501526, 511005, 525293];

m11=[197108, 257302, 270079, 272322, 277020, 283637, 284948, 290643, 292355, 294510, 296787,...
    299935, 307438, 323851, 326652, 327405, 330553, 335078, 336796, 346748, 347785, 348894,...
    349897, 350486, 350843, 363981, 369789, 392961, 394719, 398282];

m115=[0, 118679, 121863, 128393, 137885, 149169, 152031, 154493, 156809, 168390, 17185, 175238,...
    184828, 186927, 204234, 207473, 214452, 225155, 226194, 230235, 234514, 247348, 256507,...
    265267, 267323, 31342, 52618, 60731, 69507, 83280];

lim1=9:0.5:11.5;
lim2=lim1+0.5;
lim2(end)=13;


for k=1:length(lim1)
    
    fprintf('plotting in mass range 1e%s to 1e%s \n$',num2str(lim1(k)),num2str(lim2(k)))

    switch k
        case 1
            idList=m9;
            massTag='m9';
        case 2
            idList=m95;
            massTag='m95';
        case 3
            idList=m10;
            massTag='m10';
        case 4
            idList=m105;
            massTag='m105';
        case 5
            idList=m11;
            massTag='m11';
        case 6
            idList=m115;
            massTag='m115';
    end
          
    baseName=['galMap_id%i_%s_%s_' massTag '_snp%i_' simDisplayName];
    
    projTag={'zy','xz','xy'};
    
    barr='[';
    for i=1:length(idList)
        
        
        idGal=idList(i);
        
        galMass=massAllGals(idGal+1);
               
        
        idHost=subsInfo.hostFof(idGal+1);
        hostMass=fofs.Group_M_Crit200(idHost+1).*illUnits.massUnit;
        
        center=subs.SubhaloPos(:,idGal+1);
        
        rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,idGal+1);
        rhalfStars=subs.SubhaloHalfmassRadType(illustris.partTypeNum('star')+1,idGal+1);
        r200=fofs.Group_R_Crit200(idHost+1);
        r500=fofs.Group_R_Crit500(idHost+1);
        
        
        
        
        rLimitOuter=2*rhalfGas;
        rLim=1.05.*rLimitOuter;
        boxx=2*rLim;
        
        circ(1).radius=2.*rhalfStars;
        circ(1).width=1.5;
        circ(2).radius=rhalfGas;
        circ(2).width=1.5;
        circ(3).radius=r500;
        circ(3).width=1.5;
        circ(4).radius=r200;
        circ(4).width=1.5;
        
        
        
        %% load  and  plot star data
        starFields={'Coordinates','Masses','GFM_StellarFormationTime'};
        stars=illustris.snapshot.loadSubhalo(bp,snap,idGal,'stars',starFields);
        stars.newCoord=illustris.utils.centerObject(stars.Coordinates,center);
        starMask=abs(stars.newCoord(1,:))<=rLim & abs(stars.newCoord(2,:))<=rLim & abs(stars.newCoord(3,:))<=rLim ;
        starMask=starMask & stars.GFM_StellarFormationTime>=0; % remove wind particles
        
        
        %'figure',hf,'axes',axStar,'zoom',rLimitOuter,'clims',[5 9],'labels','no','arrow',scaleArw,'text',text);
        %
        
        % load gas data
        gasFields={'Coordinates','Masses','Density','ElectronAbundance'...
            'EnergyDissipation','GFM_CoolingRate','GFM_Metallicity',...
            'InternalEnergy','Machnumber','MagneticField','NeutralHydrogenAbundance',...
            'Potential','StarFormationRate','Velocities'};
        
        gas=illustris.snapshot.loadHalo(bp,snap,idHost,'gas',gasFields);
        
        gas.newCoord=illustris.utils.centerObject(gas.Coordinates,center);
        
        gas=illustris.utils.addEntropy(gas);
        gas=illustris.utils.addTemperature(gas);
        gas=illustris.utils.addPressure(gas);
        
        gasMask=abs(gas.newCoord(1,:))<=rLim & abs(gas.newCoord(2,:))<=rLim & abs(gas.newCoord(3,:))<=rLim ;
        
        %vcm=sum(gas.Velocities(:,gasMask).*gas.Masses(gasMask),2)./sum(gas.Masses(gasMask));
        
        %% plot
        ng=256;
        %titleTag=['gal ID: ' num2str(idGal)];
        
        mgalText=mk_mvir_string(galMass,'gal');
        mhostText=mk_mvir_string(hostMass,'host');
        ssfrString=mk_exponent_string(ssfr(idGal+1));
        titleTag=['gal ID: ' num2str(idGal) ', '   mgalText(1:end-1) ', ' mhostText(2:end), ', $\mathrm{sSFR}=' ssfrString '\,\mathrm{yr^{-1}}$'];
        if ssfr(idGal+1)>=1e-11
            qTag='SF';
        else
            qTag='Q';
        end
        
        %massText.str=sprintf('%s, %s',mgalText,mhostText);
        % massText.pos=[-70 70 100 10];
        % massText.color='k';
        % massText.fontsize=13;
        
        cmapStar=brewermap(256,'*Greys');
        cmapTemp=brewermap(256,'*Spectral');
        cmapDens=brewermap(256,'*YlGnBu');
        cmapEnt=brewermap(256,'PuRd');
        cmapTc=brewermap(256,'YlOrRd');
        cmapPress=brewermap(256,'YlOrBr');
        caxTemp=[4 6.5];
        caxEnt=[-1.5 1.5];
        caxDens=[-5 -1];
        caxTc=[-1 1.5];
        caxPress=[-5 -2];
        
        thickness=1.5.*rhalfGas;
        
        % run over projection
        for j=1:length(projTag)
            
            hf=figure;
            set(hf,'position',[ 25          18        1548         967],'color','k');
            
            
            %% stars
            
            axStar=axes('Parent',hf,...
                'Position',[0.09 0.53  0.26 0.35],'color','w');
            
            
            illustris.plots.mkmapStars('star',stars,'type','mass','ng',ng,'mask',starMask,'box',boxx,...
                'zoom',rLimitOuter,'clims',[5 9],'circ',circ,'noXlab',...
                'figure',hf,'axes',axStar,projTag{j});
            ht=titlemine('Stars',16);set(ht,'color','w');
            
            %% density
            %subplot(2,3,4)
            axDens=axes('Parent',hf,...
                'Position',[0.09 0.1 0.26 0.35]);
            illustris.plots.mkmapGas('gas',gas,'type','ndensity','ng',ng,...
                'mask',gasMask,'circ',circ,'clims',caxDens,...%'vfield','dilute',6,'vcm',vcm,...
                'box',boxx,'zoom',rLimitOuter,'thick',thickness,'nobackground',...
                'figure',hf,'axes',axDens,projTag{j},'black');  %,...%'streamdense',0,...
            
            ht=titlemine('Density',16);set(ht,'color','w');
            
            
            %% entropy
            axEnt=axes('Parent',hf,...
                'Position',[0.4 0.53  0.26 0.35]);
            
            illustris.plots.mkmapGas('gas',gas,'type','entropy','ng',ng,...
                'mask',gasMask,'circ',circ,'clims',caxEnt,'noXlab','noYlab',...%'vfield','dilute',6,'vcm',vcm,...
                'box',boxx,'zoom',rLimitOuter,'thick',thickness,'nobackground',...
                'figure',hf,'axes',axEnt,projTag{j},'black');
            
            ht=titlemine('Entropy',16);set(ht,'color','w');
            
            %% temp
            %%subplot(2,3,2)
            axTemp=axes('Parent',hf,...
                'Position',[0.4 0.1  0.26 0.35]);
            
            illustris.plots.mkmapGas('gas',gas,'type','temp','ng',ng,...
                'mask',gasMask,'circ',circ,'clims',caxTemp,'noYlab',...%'vfield','dilute',6,'vcm',vcm,...
                'box',boxx,'zoom',rLimitOuter,'thick',thickness,'nobackground',...
                'figure',hf,'axes',axTemp,projTag{j},'black'); %,...%'streamdense',0,...
            
            ht=titlemine('Temperature',16);set(ht,'color','w');
            
            
            %% Pressure
            axPres=axes('Parent',hf,...
                'Position',[0.7 0.53  0.26 0.35]);
            
            illustris.plots.mkmapGas('gas',gas,'type','pressure','ng',ng,...
                'mask',gasMask,'circ',circ,'noXlab','noYlab','clims',caxPress,...%'vfield','dilute',6,'vcm',vcm,...
                'box',boxx,'zoom',rLimitOuter,'thick',thickness,'nobackground',...
                'figure',hf,'axes',axPres,projTag{j},'black'); %,...%'streamdense',0,...
            
            ht=titlemine('Pressure',16);set(ht,'color','w');
            
            %% tcool
            axTcool=axes('Parent',hf,...
                'Position',[0.7 0.1 0.26 0.35]);
            illustris.plots.mkmapGas('gas',gas,'type','tcool','ng',ng,...
                'mask',gasMask,'circ',circ,'clims',caxTc,'noYlab',...%'vfield','dilute',6,'vcm',vcm,...
                'box',boxx,'zoom',rLimitOuter,'thick',thickness,'zy','nobackground',...
                'figure',hf,'axes',axTcool,projTag{j},'black'); %,...%'streamdense',0,...
            
            ht=titlemine('Cooling Time',16);set(ht,'color','w');
            
            %% set colors
            
            colormap(axStar,cmapStar);
            colormap(axTemp,cmapTemp);
            colormap(axDens,cmapDens);
            colormap(axEnt,cmapEnt);
            colormap(axTcool,cmapTc);
            colormap(axPres,cmapPress);
            
            
            zstr=titleTag;  %sprintf('$z=%1.2f $',zred(indx));
            annotation(hf,'textbox',...
                [0.18 0.95 0.8 0.05],...
                'string',zstr,'Fontsize',20,...
                'Interpreter','latex',...
                'color',[1 1 1],'linestyle','none')
            
            
            fname=sprintf(baseName,idGal,qTag,projTag{j},snap);
            printout_fig(hf,fname,'subdir','mapSurvey','nopdf')
            
            close(hf)
            
        end
        
        barr=[barr '='];
        
        if i==10 || i==20
            barr=[barr '|'];
        elseif i==30
             barr=[barr ']\n'];
        end
        
        fprintf(barr);
    end
    
    
end

%
%
%
% hf=illustris.plots.mkmapGas('gas',gas,'type','entropy','ng',ng,...
%     'mask',gasMask,'circ',circ,'brewer',brewEnt,'clims',caxEnt,...%'vfield','dilute',6,'vcm',vcm,...
%     'box',boxx,'zoom',rLimitOuter,'thick',rhalfGas,'nobackground',...
%     'text',massText,'title',titleTag); %,...%'streamdense',0,...
% for j=1:length(hf)
%     fname=sprintf(baseName,'ent',idGal,projTag{j},snap);
%     printout_fig(hf(j),fname,'subdir','mapSurvey','v','png')
% end



%'labels','no','black','figure',hf)

%printout_fig(gcf,'testTemp')
%set(gca,'XTickLabel','','YTickLabel','')

%ht=titlemine('Gas Temperature',20);set(ht,'color','w');
