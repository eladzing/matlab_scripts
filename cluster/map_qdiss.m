

%tag1='test_dis_sm150';
%printoutdir='/home/eladzing/Ubuntu One/cluster/printout/dissipation_tests';
printoutdir='C:/Users/owner/Ubuntu One/cluster/printout/dissipation_tests';

for kk=1:2
    switch kk
        case 1
            cl=[-10 -2];
            unitag='$\log E\,[\mathrm{\frac{G M_{vir}}{R_{vir}}}]$';
            tag1='testdis_efull';
        case 2
            cl=[32 38];
            unitag='$\log \frac{E}{m}\,[\mathrm{\frac{G}{R_{vir}}}]$';
            tag1='testdis_espec';
    end
    for k=1:6
        
        switch k
            case 1
                dcub=qdiss(:,:,:,kk);tag='qdis';
            case 2
                dcub=deltaK(:,:,:,kk);tag='dkin';
            case 3
                dcub=deltaP(:,:,:,kk);tag='dpot';
            case 4
                dcub=deltaT(:,:,:,kk);tag='dth';
            case 5
                dcub=pdv(:,:,:,kk);tag='pdv';
            case 6
                dcub=qcube(:,:,:,kk);tag='cool';
        end
        pmask=dcub>0;
        
        mkmap('data',dcub.*pmask,'clims',cl,'proj',[1 0 0],'normalize',normfac(kk),...
            'title',sprintf('$Q_{%s}>0$',tag),'titletype','custom','boxx',boxx,'log','on','marks',[0.1 0.2 0.5]*get_rvir,'bartag',unitag,...
            'print','yes','format','both','printag',sprintf('%s_%s_pos',tag,tag1),'printout',printoutdir);
        if(~strcmp(tag,'cool'))
            mkmap('data',-1.*dcub.*~pmask,'clims',cl,'proj',[1 0 0],'normalize',normfac(kk),...
                'title',sprintf('$Q_{%s}<0$',tag),'titletype','custom','boxx',boxx,'log','on','marks',[0.1 0.2 0.5]*get_rvir,'bartag',unitag,...
                'print','yes','format','both','printag',sprintf('%s_%s_neg',tag,tag1),'printout',printoutdir);
        end
    end
    %pause
    %close all
    cl=[1 1];
    
    switch kk
        case 1
            et='E';
        case 2
            et='e';
    end
    
    
    mkmap('data',squeeze(ekin(:,:,:,kk)),'clims',cl,'proj',[1 0 0],'normalize',normfac(kk)','title',sprintf('$ %s_{k}$',et),'boxx',boxx,'log','on','bartag',unitag,...
        'print','yes','format','both','printag',sprintf('%s_%s','kin',tag1),'printout',printoutdir);
    mkmap('data',abs(squeeze(egrav(:,:,:,kk))),'clims',cl,'proj',[1 0 0],'normalize',normfac(kk)','title',sprintf('$ %s_{p}$',et),'boxx',boxx,'log','on','bartag',unitag,...
        'print','yes','format','both','printag',sprintf('%s_%s','pot',tag1),'printout',printoutdir);
    mkmap('data',squeeze(eth(:,:,:,kk)),'clims',cl,'proj',[1 0 0],'normalize',normfac(kk)','title',sprintf('$ %s_{th}$',et),'boxx',boxx,'log','on','bartag',unitag,...
        'print','yes','format','both','printag',sprintf('%s_%s','eth',tag1),'printout',printoutdir);
end
close all
%mkmap('data',abs(qdiss),'clims',[-6 6],'normalize',normfac,'proj',[1 1 1],'title','$Q_{dis}$','boxx',boxx,'log','off','marks',[0.1 0.2 0.5]*get_rvir);
%mkmap('data',abs(delitaK),'clims',[-6 6],'normalize',normfac,'proj',[1 1 1],'title','$Q_{k}$','boxx',boxx,'log','off','marks',[0.1 0.2 0.5]*get_rvir);
%mkmap('data',abs(deltaP),'clims',[-6 6],'normalize',normfac,'proj',[1 1 1],'title','$Q_{p}$','boxx',boxx,'log','off','marks',[0.1 0.2 0.5]*get_rvir);
%mkmap('data',abs(qcube),'clims',[-6 6],'normalize',normfac,'proj',[1 1 1],'title','$Q_{cool}$','boxx',boxx,'log','off','marks',[0.1 0.2 0.5]*get_rvir);
%mkmap('data',abs(pdv),'clims',[-6 6],'normalize',normfac,'proj',[1 1 1],'title','$Q_{pdv}$','boxx',boxx,'log','off','marks',[0.1 0.2 0.5]*get_rvir);