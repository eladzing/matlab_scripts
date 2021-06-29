list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%list1=[101 24];

result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout/';
set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);

pflag='print';
hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    

rmf=[0 0.02 0.1 0.6];
rxf=[0.01 0.1 0.3 1.2];

ptags={'core_clean','inner_clean','inner_clean2','outer_clean'};


len=1000;

smin=-2;smax=2;
for ii=1:4
  rmfac=rmf(ii);
  rxfac=rxf(ii);
  printag=char(ptags(ii));   
  
  
for id=1:length(list1)
    halopath=sprintf(hpath,list1(id));
    clustername=sprintf('CL%d',list1(id))
    stack{id,1}=clustername;
    switch list1(id)
        case {101,102,103,104,105,106,107,5}
            smallbox=2;bigbox=8;
        otherwise
            smallbox=1;bigbox=8;
    end
    smass_hist 
    stack{id,1}=clustername;
    stack{id,2}=tax;
    stack{id,3}=sh_in(:,1)/dt;
    stack{id,4}=sh_out(:,1)/dt;
 end
      
  save(sprintf('mat_files/stk_shist_%s.mat',printag),'stack');
  
end  
  

%[stxin stdin]=stax(stack,tax,2,3,9);
%[stxou stdou]=stax(stack,tax,2,4,9); 

%figure;plot(tax,[stxin' stxou']); grid; xlim([mnt mxt]);
%set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
%xlabel('log T/T_{vir}');
%ylabel('d(M/Mgas)/dlog(T/T_{vir})');
%title(sprintf( '%s Infolw/Outflow Temperature Histogram (%s<r/R_{vir}<%s)','Stacked',num2str(rmfac,2),num2str(rxfac,2)));
%legend('Inflow', 'Outflow','Location','NorthWest');

%saveas(gcf,sprintf('%s/%s_thist_%s.png',result_dir,'stack',printag));


