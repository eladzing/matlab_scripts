q10=gradesallunified2.Q106>0;
q11=gradesallunified2.Q116>0;
q17=gradesallunified2.Q1710>0;
q18=gradesallunified2.Q1810>0;
qCountTot=q10+q11+q17+q18;

%% 

qCount=zeros(size(IDNumber));
for j=1:length(IDNumber)
    ind=find(IDNumber(j)==gradesallunified2.IDNumber);
    if length(ind)>1
        error('Something weird %i',j)
    end
    qCount(j)=qCountTot(ind);
end

%% 
fid=fopen('C:\Users\eladz\Documents\teach\physics_for_bios_77133\moed_A\question_count_list_full.csv','w');
for i=1:length(IDNumber)
    fprintf(fid,'%i8 , %i, \n',IDNumber(i),qCount(i));
end
fclose(fid);

fid=fopen('C:\Users\eladz\Documents\teach\physics_for_bios_77133\moed_A\question_count_list.txt','w');
for i=1:length(IDNumber)
    st=num2str(IDNumber(i));
    fprintf(fid,'%s , %i, \n',st(end-3:end),qCount(i));
    id2=str2num(st(end-3:end));
end
fclose(fid);


%% test wit;h reuven 

for i=1:length(IDNumber)
    st=num2str(IDNumber(i));
    id2(i)=str2double(st(end-3:end));
end


rev=[0034, 	3;0042, 	3;0056, 	3;
0066,	2;0078, 	3;0100,	4;0249,	1;0271, 	2;0284, 	3;0314, 	3;0323, 	0;0352, 	2;0594, 	3;0615, 	3;0711, 	2;0796, 	3;0804, 	2;0846, 	3;1068, 	2;1103, 	3;1171, 	2;1220, 	4;1241, 	3;1290, 	4;1315, 	3;1324, 	3;...
1423, 	3;1464, 	2;1606, 	2;1998, 	2;2054, 	1;2406, 	0;2434, 	3;2650, 	3;2658, 	2;2674, 	4;2888, 	4;2901, 	2;2925, 	4;2926, 	4;2930, 	1;3114, 	2;3126, 	1;3182, 	2;...
3183, 	3;3247, 	3;3263, 	3;3308, 	1;3353, 	2;3355, 	3;3360, 	3;3414, 	1;3512, 	2;3585, 	3;3599, 	3;3636, 	2;3652, 	4;3681, 	4;3817, 	4;3820, 	4;3863, 	3;3896, 	4;...
3979, 	2;4031, 	2;4153, 	3;4217, 	2;4227, 	3;4303, 	2;4306, 	3;4609, 	2;4651, 	3;4766, 	3;4850, 	1;4852, 	2;4898, 	4;4945, 	1;4992, 	1;5136, 	1;5201, 	2;...
5449, 	1;5459, 	3;5558, 	2;5826, 	2;5876, 	1;5931, 	1;5961, 	3;6128, 	2;6153, 	2;6166, 	3;6299, 	4;6326, 	3;6346, 	2;6488, 	3;6647, 	4;6669, 	3;6676, 	3;6719, 	4;...
6723, 	2;6881, 	3;6926, 	3;6971, 	3;6975, 	1;7031, 	3;7198, 	3;7452, 	3;7465, 	2;7555, 	3;7557, 	2;7587, 	2;7724, 	3;...
7736, 	4;7750, 	4;7817, 	2;7830, 	2;7976, 	3;8016, 	3;8093, 	2;8192, 	1;8306, 	3;8487, 	4;8539, 	2;8573, 	4;...
8703, 	2;8734, 	2;8794, 	3;8928, 	4;8999, 	3;9002, 	2;9033, 	3;9119, 	3;9202, 	2;9243, 	3;9268, 	4;9298, 	2;...
9342, 	3;9385, 	3;9479, 	2;9513, 	3;9608, 	3;9650, 	1;9706, 	2;9726, 	3;9792, 	2;9835, 	3;9916, 	2;9943, 	0;9954,	2];


for i=1:length(id2)
    ii=find(id2(i)==rev(:,1));
    test(i)=qCount(i)==rev(ii,2);
end