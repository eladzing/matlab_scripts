list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
r=0.01:0.01:10;
figure
for i=1:length(list)
    new_env(sprintf('CL%d',list(i)),'csf','a1','win')

    s=read_S_Profile(r);
    
    if i>1
        hold on
    end
    loglog(r/0.7*1000,s)
end 
xlabelmine('$r\,[\mathrm{kpc}]$')
ylabelmine('$K\,[\mathrm{KeV\,cm^2}]$')