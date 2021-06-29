%% try out Andrey's suggestion for counting the total angle of inflow
% above a given limit
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];


 dsProf=zeros(639,length(list));
        rrProf=dsProf;
for k=1:length(list)
    new_env(list(k));
    
    %if k==1
        global NCELL
        global hub
        
       
        
        r0=(1:NCELL)-0.5;
        rrSph0=zeros(NCELL,NCELL,NCELL);
        
        for i=1:NCELL
            rrSph0(i,:,:)=r0(i);
        end
    %end
    rr=zeros(NCELL,4);
    dsInflow=zeros(NCELL,4);
    for i=1:4
        boxx=floor(2.^(i-1));
        fluxnorm=(0.1117.*(get_mvir./1e15).^0.15*(1+0)^2.25)./(4.*pi); % normalization for flux
        %% load cube
        
        %% calculate inflow on sphere
        
        fluxShell=flux_sphere(boxx)./fluxnorm;
        inflowMask=fluxShell<=-5;
        
       % if k==1
            rrSph=rrSph0.*(0.5*boxx/hub/NCELL);
            rr(:,i)=r0.*(0.5*boxx/hub/NCELL);
            
        %end
        dOmega=ds_sphere(boxx)./rrSph.^2;
        
        dsInflow(:,i)=sum(sum(dOmega.*inflowMask,3),2)./(4*pi);
        
    end
    
    dsProf(:,k)=cat(1,dsInflow(1:255,1),dsInflow(128:255,2),dsInflow(128:255,3),dsInflow(128:255,4));
    rrProf(:,k)=cat(1,rr(1:255,1),rr(128:255,2),rr(128:255,3),rr(128:255,4))./get_rvir;

%     hf=figure;
%     semilogx(rr(:,1),dsInflow(:,1),rr(:,2),dsInflow(:,2),rr(:,3),dsInflow(:,3),rr(:,4),dsInflow(:,4))
%     hold on
%     semilogx(rrProf(:,k),dsProf(:,k))
%     pause
%     close(hf) 
end


figure

for i=1:16
    h(i)=semilogx(rrProf(:,i),dsProf(:,i),'DisplayName',num2str(list(i)));
    %if i==1
        hold on
    %end
    
end
xlim([0.05 1]);
ylim([0 1])

xlabelmine('$r/R_{\mathrm{vir}}$');
ylabelmine('$\mathrm{d}\Omega_{\mathrm{inflow}}/4\pi$')

clickableLegend(h)
%hl=legend(h);
%set(hl,'Interpreter','latex','Fontsize',14)


% %% find density histogram
%     roShell=RHOG_sphere(boxx);
%
%
%
%
%     if boxx==1
%        startShell=1;
%     else
%     startShell=129;
%     end
%
%         for j=startShell:size(fluxShell,1)-1
%
%         end
% end

%figure
%semilogx(rr(1:255,1),dsInflow(1:255,1),rr(128:255,2),dsInflow(128:255,2),rr(128:255,3),dsInflow(128:255,3),rr(128:255,4),dsInflow(128:255,4))