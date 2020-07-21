list=[103 106 6 11 14 ];
%list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
ind=[ 3  6  10  14 15];
%ind=1:length(list);
prj={'xz' 'xy' 'all' 'xz' 'xy'};
load('/home/zinger/workProjects/matlab_scripts/cluster/mat_files/shockedge.mat')


units;
%prj='all';
for i=1 %:length(list)
    
    if list(i)==6
        labs='full';
    else
        labs='half';
    end
    
    
    for j=1:2
        switch j
            case 1
                a='a1';
                shk=shockedge_a1;
                shkM1=shockedgeMask1_a1;
                shkM2=shockedgeMask2_a1;
            case 2
                a='a06';
                shk=shockedge_a06;
                shkM1=shockedgeMask1_a06;
                shkM2=shockedgeMask2_a06;
        end
        
        new_env(list(i),a)
        
        c1=shk{ind(i),4};
        c2=shk{ind(i),5};
        
        circ=cat(2,c1(shkM1(ind(i),:)),c2(shkM2(ind(i),:)));
        circ2=cat(2,c1,c2);
        
        alf=0.5;
        ro=RHOG(8);
        rps=log10(alf.*(ro.*(Units.Ms/Units.Mpc^3)).*(get_vvir.*Units.km).^2);
        
        mkmap(8,'data',rps,'labels',labs,'brewer','YlOrRd',prj{i},'vfield',...
            'vfieldcolor','k','streamlines','off','dilute',10,'bartag','$\log\,P_{\mathrm{ram}}\,[\mathrm{erg\,cm^{-3}}]$',...
            'title','none','marks',circ,'proper','clims',[-15 -10.5],'print','printag','papIIRev2');
        % %     mkmap(8,'type','s','clims',[0 5],'labels',labs,'brewer','*RdYlBu',prj,...
        %          'title','none','marks',circ2,'noprint','printag','papII','proper');
        %     mkmap(8,'type','t','clims',[4 8],'labels',labs,'brewer','*RdBu',prj{i},...
        %         'title','none','marks',circ,'print','printag','papIIRev','proper');
        %     mkmap(8,'type','n','labels',labs,'clims',[-7 -2],'brewer','*RdBu',prj{i},...
        %         'title','none','marks',circ,'print','printag','papIIRev','proper');
        %   %'clims',[4 8],
        
        %     mkmap(8,'type','pressure','labels',labs,'brewer','*RdYlBu',prj{i},...
        %         'title','none','marks',circ,'print','printag','papII','proper');
        
        %close all;
    end
end
