for i=iii'
    
    figure(4)
    tim=redshift2time(outskirtJF_histories.hist(i).zr,'cosmo',cosmoStruct).lookback;
    plot(tim,outskirtJF_histories.hist(i).rpos);
    hold on
    plot(tim,outskirtJF_histories.hist(i).isSat);
    hold off
    titlemine(num2str(i));
   
    pause
    
end

%%  calculate dynamical time 
rhoc=rho_crit(zred,'cosmo',cosmoStruct);

ts=sqrt((3/(4*pi*5/3)).*(Units.G.*200*rhoc'.*Units.Ms./Units.Mpc^3).^-1)./Units.Gyr;
tdyn=sqrt((3*pi/32).*(Units.G.*200*rhoc'.*Units.Ms./Units.Mpc^3).^-1)./Units.Gyr;


%% 
spb_minTime=10;
spb_inTime=1.5;
spb_inRad=10;

spbNum=sum(outskirtJF_histories.rposMin<spb_inRad & ...
    outskirtJF_histories.timeMin<spb_minTime & ...
    outskirtJF_histories.timeLastInRv<spb_inTime )
    
%%

rr=0.1:0.05:1;

for i=1:length(outskirtJF_histories.hist)
    zr=outskirtJF_histories.hist(i).zr;
    for j=1:length(rr)
        zlast=zr(find(outskirtJF_histories.hist(i).rpos(outskirtJF_histories.hist(i).isSat)<rr(j),1,'first'));   
        if ~isempty(zlast)
        outskirtJF_histories.lastTime(i,j)=redshift2time(zr(1),'cosmo',cosmoStruct).age-...
            redshift2time(zlast,'cosmo',cosmoStruct).age;  % time since last time it was within r200 of host.
        else
            outskirtJF_histories.lastTime(i,j)=-1;
        end
            
    end
end

%% normalize by sound crossing time 
lastt=outskirtJF_histories.lastTime;
for i=1:length(rr)
    lastt(:,i)=outskirtJF_histories.lastTime(:,i)./ts;
end
lastt(lastt<0)=10;

%% 
tt=1:0.25:3;

spb_num=zeros(length(tt),length(rr));

for i=1:length(tt)
    for j=1:length(rr)
        
        spb_num(i,j)=sum(lastt(:,j)<=tt(i))./1168;
    end
end


    


%%
spb_inTime=1.5;
spbNum=sum(outskirtJF_histories.lastTime(i,1)<spb_inRad & ...
    outskirtJF_histories.timeMin<spb_minTime & ...
    outskirtJF_histories.timeLastInRv<spb_inTime )
    
