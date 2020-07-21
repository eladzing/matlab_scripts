aexp='a1';typ='csf';
hubbleflag='hub';
dilute = 8;
clims=[-3 3];
%result_dir='/home/eladzing/work/sshfs/titan3/cold_flows/tit3/printout/';
result_dir='/home/eladzing/work/cold_flows/datacube/printout/';
printag='zonepen_new';
pflag='print';

list=[3 11 14 6];

%% CL3
for id=1:length(list)
  clustername=sprintf('CL%d',list(id));
  new_env(clustername,typ,aexp);
  switch id
      case {1, 2, 3}
         plotproj=[0 0 1];
         prjtag='XY'
      case 4
         plotproj=[1 0 0];
         prjtag='YZ'
  end
         
  for bx=0:3
    
    boxx=2.^bx
    %clims=[-3 3];
    figure;

    zone_pen_block(plotproj,boxx,clims,hubbleflag);

    if strcmp(pflag,'print')
      saveas(gcf,sprintf('%s/%s_map%s_b%d_%s.png',result_dir,clustername,prjtag,boxx,printag)); 
      %set(gcf, 'PaperPositionMode', 'auto');
      
      print('-depsc2',sprintf('%s/%s_map%s_b%d_%s.eps',result_dir,clustername,prjtag,boxx,printag));
      %saveas(gcf,sprintf('%s/%s_map%s_b%d_%s.eps',result_dir,clustername,prjtag,boxx,printag));    
    end

  end
end
%clear all
%close all


%% old




% boxx=1;
% figure;
% clims=[-3 3];
% 
% zone_pen_block(plotproj,boxx,clims,hubbleflag);
% 
% if strcmp(pflag,'print')
%     saveas(gcf,sprintf('%s/%s_map%s_b%d_%s.png',result_dir,clustername,prjtag,boxx,printag));    
% end
%clear all
% 
% %======================
% %% CL14
% zone_pen_param
% clustername='CL14';
% new_env(clustername,typ,aexp);
% plotproj=[0 0 1];
% 
% my_subplot(4,2,3);
% boxx=8;
% %clims=[-3 3];
% 
% zone_pen_block;
% ylabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
% %----------------------
% 
% my_subplot(4,2,4);
% boxx=1;
% %clims=[-3 3];
% 
% zone_pen_block;
% ylabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
% 
% bar=colorbar;
% set(get(bar,'Title'),'String','$\frac{dot M/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$','Fontsize',12,'Interpreter','latex');
% 
% clear all
% %======================
% %% CL6
% zone_pen_param
% clustername='CL6';
% new_env(clustername,typ,aexp);
% plotproj=[1 0 0];
% 
% my_subplot(4,2,5);
% boxx=8;
% %clims=[-3 3];
% 
% zone_pen_block;
% ylabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
% %----------------------
% 
% my_subplot(4,2,6);
% boxx=1;
% %clims=[-3 3];
% 
% zone_pen_block;
% 
% bar=colorbar;
% set(get(bar,'Title'),'String','$\frac{dot M/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$','Fontsize',12,'Interpreter','latex');
% 
% clear all
% %======================
% %% CL11
% zone_pen_param
% clustername='CL11';
% new_env(clustername,typ,aexp);
% plotproj=[0 0 1];
% 
% my_subplot(4,2,7);
% boxx=8;
% %clims=[-3 3];
% 
% zone_pen_block;
% xlabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
% ylabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
% %----------------------
% 
% my_subplot(4,2,8);
% boxx=1;
% %clims=[-3 3];
% 
% zone_pen_block;
% 
% bar=colorbar;
% set(get(bar,'Title'),'String','$\frac{dot M/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\left|_{\mathrm{virial}}]$','Fontsize',12,'Interpreter','latex');
% 
% xlabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
% 
% %======================
% 
% if strcmp(pflag,'print')
%   saveas(gcf,sprintf('%s/flux_map_multi_pap.png',result_dir));    
%   saveas(gcf,sprintf('%s/flux_map_multi_pap.eps',result_dir));
% end
% 
% clear all
