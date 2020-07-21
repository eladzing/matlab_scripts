list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
typ='csf';% 'adiabatic');
aexp='a1';

%cvir0=10;

ngrid=256;
max=4./0.7;
min=1./0.7./ngrid;
np=max/min-1;

rmin=[min log10(min)];
rmax=[max log10(max)];
dr=2.*([min ((rmax(2)-rmin(2))./np)]);

rp=rmin(1):dr(1):rmax(1);
lrp=rmin(2):dr(2):rmax(2); 
lrp=10.^lrp;

rax=zeros(2,size(rp,2));
rhog=zeros(2,size(rp,2));
rhodm=zeros(2,size(rp,2));
rhot=zeros(2,size(rp,2));

rax(1,:)=rp;
rax(2,:)=lrp;

%lrpn=[];
%lrot=[];
%lrod=[];

for id=1:length(list1)
    clustername=sprintf('CL%d',list1(id));
    new_env(clustername,typ,aexp); 
    rho_profs{id,1}=clustername;
    %global hub 
    
    rv=get_rvir();
    mv=get_mvir();  
    
    %rmin=rv.*0.008;
    %rmax=rv.*2.0;
    %dr=(rmax-rmin)./1000;
    
    %rp=rmin(1):dr(1):rmax(1);
    %lrp=rmin(2):dr(2):rmax(2); 
    %lrp=10.^lrp;
    %rpl=-2:0.1:log10(4/hub);
    %rpl=10.^rpl;
    
    [rog, rot] =read_RHO_Profiles(rp);
    [lrog, lrot] =read_RHO_Profiles(lrp);  
    mt=read_MTOT_Profile(rp);
    lmt=read_MTOT_Profile(lrp);
    
    rhog(1,:)=rog;
    rhog(2,:)=lrog;
    rhot(1,:)=rot;
    rhot(2,:)=lrot;
    rhodm(1,:)=rot-rog;
    rhodm(2,:)=lrot-lrog;
    mtot(1,:)=mt;
    mtot(2,:)=lmt;
    
    rho_profs{id,2}=rax; %!./rv;
    rho_profs{id,3}=rhot; 
    rho_profs{id,4}=rhodm; 
    rho_profs{id,5}=rhog; 
    rho_profs{id,6}=mtot;
    rho_profs{id,7}=[rv,mv];
    
    
    %rho_profs{id,2}=rp; %!./rv;
    %rho_profs{id,3}=lrp; %!./rv;
    %rho_profs{id,4}=rot;
    %rho_profs{id,5}=lrot;
    %rho_profs{id,6}=rog;
    %rho_profs{id,7}=lrog;
    %rho_profs{id,8}=rot-rog;%rodm;
    %rho_profs{id,9}=lrot-lrog;%rodm;
    %rho_profs{id,10}=[rv,mv];
    %lrods=lrod-1;
    
    %ind=find(rpn>=1/cvir0 & rpn<1.0);
    
    %roz=4*rott(ind(1));
    
      
    
    %clname{id}=sprintf('CL%s',num2str(flist(id)));
    
    %ft = fittype( 'rho_nfw(r,ro0,cvir)');
    %f = fit(rpn,rott,ft, 'StartPoint', [roz, cvir0]);
    
    
    %rot(end+1,:)=rott;
    %rpn(end+1,:)=rp./rv;
    
    
end
save(sprintf('mat_files/rho_profs_%s.mat',aexp),'rho_profs');
%plot( f, x, y ) 


%figure('PaperType','A4','PaperOrientation','landscape','PaperPosition',[0.5 0.5 10.5 8]);
    %title(sprintf('%s Clusters Mass Flux',titletag),'Fontsize',12);
%     for j=1:length(list1)
%         len=0.22;
%         a=mod(j-1,4);
%         b=4-ceil(j/4);
%         left=0.1+len.*a ;
%         bottom= 0.1+len.*b;
%         subplot('Position',[left bottom (len-0.02) (len-0.02)])
%         %subplot(4,4,j);
%         loglog(rpn(j,:),rot(j,:),rpn(j,:),ronfw2(j,:),'linewidth',1);grid;xlim([1e-2 6]);%;ylim([-0.38 0.2]);
%         %set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'Fontsize',12,'linewidth',1);
%         line([1e-3 1e1],[0 0],'Color','k');
%         %text(0.021,0.11,sprintf('%s',clname{j}),'Fontsize',9);
%         switch (a)
%             case 0
%                 %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',10,'Interpreter','latex');
%                 ylabel('rho','Fontsize',10,'Interpreter','latex');
%                 set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'Fontsize',12,'linewidth',1,...
%                     'xticklabel','');
%             otherwise
%                 set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'Fontsize',12,'linewidth',1,...
%                     'yticklabel','','xticklabel','');
%         end
%         switch (b)
%             case 0
%                 xlabel('$log(r/R_{\mathrm{vir}})$','Fontsize',10,'Interpreter','latex');
%                 %ylabel('$\dot{M}/M_{g} [\mathrm{Gyr}^{-1}]$','Fontsize',10,'Interpreter','latex');
%                 set(gca,'xticklabel',[0 -2 -1]);
%             otherwise
%                 %set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'Fontsize',12,'linewidth',1,...
%                 %    'yticklabel','','xticklabel','');
%         end
%     end
%     
    