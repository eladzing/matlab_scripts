list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
mv1=[];
for i=1:16
new_env(sprintf('CL%d',list(i)),'csf','a1','win')
mv1(end+1)=get_mvir;
end