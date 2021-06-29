%% penenetration vs. entropy 

penetration_depth_a1
load('C:\Users\eladzing\OneDrive\cluster\matlab\mat_files\profiles.mat')


%% entropy vs. penetration 
%penID=find(maxpen<=median(maxpen));
rpP=[];
spP=[];
rpN=[];
spN=[];

rp=0.01:0.01:4;

for i=1:16
    
    new_env(profile(1).cluster(i).name);
    if penet(i)
    
        rpP(end+1,:)=rp;
        spP(end+1,:)=read_S_Profile(rp.*get_rvir);
        %spP(end+1,:)=profile(1).cluster(i).sProf;
    else
        rpN(end+1,:)=rp;
        spN(end+1,:)=read_S_Profile(rp.*get_rvir);
        %spN(end+1,:)=profile(1).cluster(i).sProf;
    end
end

msP=mean(spP,1);
msN=mean(spN,1);
stN=std(spN,1);
stP=std(spP,1);

figure(1)
hold off

h=loglog(rp,msP,'-r',rp,msN,'-b','linewidth',2);
hold on

%[hP,m]=shade2curves(rp,msP-0.5.*stP,msP+0.5.*stN,'b','k',0.1);
%[hN,m]=shade2curves(rp,msN-0.5.*stN,msN+0.5.*stN,'r','k',0.1);
loglog(rp,msP-0.5.*stP,'--r',rp,msN-0.5.*stN,'--b')
loglog(rp,msP+0.5.*stP,'--r',rp,msN+0.5.*stN,'--b')
xlim([0.05 1])
ylim([100 3e4])
set(gca,'fontsize',12)

xlabelmine('$r/R_{vir}$')
ylabelmine('$ K \,[\mathrm{keV\,cm^2}]$');
hl=legend(h,'Deep Penetrating','No Deep Penetration');
set(hl,'Interpreter','latex','Fontsize',12,'Location','NorthWest' );

%% entropy vs. rlx 
%penID=find(maxpen<=median(maxpen));
rpR=[];
spR=[];
rpU=[];
spU=[];

rp=0.01:0.01:4;

for i=1:16
    
    new_env(profile(1).cluster(i).name);
    if rlx(i)
    
        rpR(end+1,:)=rp;
        spR(end+1,:)=read_S_Profile(rp.*get_rvir);
        %spP(end+1,:)=profile(1).cluster(i).sProf;
    else
        rpU(end+1,:)=rp;
        spU(end+1,:)=read_S_Profile(rp.*get_rvir);
        %spN(end+1,:)=profile(1).cluster(i).sProf;
    end
end

msR=mean(spR,1);
msU=mean(spU,1);
stU=std(spU,1);
stR=std(spR,1);

figure(2)
hold off

h=loglog(rp,msU,'-r',rp,msR,'-b','linewidth',2);
hold on

%[hP,m]=shade2curves(rp,msP-0.5.*stP,msP+0.5.*stN,'b','k',0.1);
%[hN,m]=shade2curves(rp,msN-0.5.*stN,msN+0.5.*stN,'r','k',0.1);
loglog(rp,msU-0.5.*stU,'--r',rp,msR-0.5.*stR,'--b')
loglog(rp,msU+0.5.*stU,'--r',rp,msR+0.5.*stR,'--b')
xlim([0.05 1])
ylim([100 3e4])
set(gca,'fontsize',12)

xlabelmine('$r/R_{vir}$')
ylabelmine('$ K \,[\mathrm{keV\,cm^2}]$');
hl=legend(h,'Unrelaxed','Relaxed');
set(hl,'Interpreter','latex','Fontsize',12,'Location','NorthWest' );

