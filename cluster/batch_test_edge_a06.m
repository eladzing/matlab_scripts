list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%ind=[ 3  6  10  14 15];
prj='all' ;%{'xz' 'xy' 'all' 'xz' 'xy'};
load('C:\Users\eladzing\Documents\cluster\matlab\mat_files\shockedge.mat')


%for i=1:length(list)
%i=1;

a='a06';
shk=shockedge_a06;

new_env(list(i),a)

c1=shk{i,4};
c2=shk{i,5};

circ=cat(2,c1(shockedgeMask1_a06(i,:)),c2(shockedgeMask2_a06(i,:)));
%circ=cat(2,c1,c2);
mkmap(8,'type','s','clims',[0 5],'labels','full','brewer','*RdYlBu',prj,...
    'marks',circ,'proper'); %'print','printag','papII','proper');
%
%     mkmap(8,'type','t','clims',[4 8],'labels',labs,'brewer','*RdYlBu',prj{i},...
%         'title','none','marks',circ,'print','printag','papII','proper');
%
%     mkmap(8,'type','pressure','labels',labs,'brewer','*RdYlBu',prj{i},...
%         'title','none','marks',circ,'print','printag','papII','proper');

%end

