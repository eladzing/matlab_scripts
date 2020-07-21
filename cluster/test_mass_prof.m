%function  test_prof(cluster)

result_dir='/home/eladzing/work/cold_flows/datacube/printout';

cluster='CL101';
type='csf';
aexp='a1';

new_env(cluster,type,aexp);

rp=0:0.001:10;

[rog rot] =read_RHO_Profiles(rp);
mt = read_MTOT_Profile(rp);

rott=rot(rot>0 & mt>0);
rpr=rp(rot>0 & mt>0);
mtt=mt(rot>0 & mt>0);




mt2=[];
r2=[];
mt3=[];
for i=2:length(rpr)
    r2(i-1)=(rpr(i)+rpr(i-1))*0.5;
    switch i
        case 2
            mt2(i-1)=4.*pi.*rott(i-1).*r2(i-1).^2.*(rpr(i)-rpr(i-1));
        otherwise
            mt2(i-1)=4.*pi.*rott(i-1).*r2(i-1).^2.*(rpr(i)-rpr(i-1))+mt2(i-2);
    end
end

for i=1:length(rpr)
    %r2(i-1)=(rpr(i)+rpr(i-1)).*0.5;
    switch i
        case 1
            mt3(i)=4.*pi.*rott(i).*rpr(i).^3;
        otherwise
            mt3(i)=4.*pi.*rott(i).*rpr(i).^2.*(rpr(i)-rpr(i-1))+mt3(i-1);
    end
end
% 
  %figure;
  %loglog(rpr,mtt,'b',r2,mt2,'r',rpr,mt3,'g');
  
  rro=rpr.^2.*rott.*4.*pi;
  mt4=cumtrapz(rpr,rro);
  
  r3=rpr(1):0.0001:10;

  m1=spline(rpr,mtt,r3);
  m2=spline(r2,mt2,r3);
  m3=spline(rpr,mt3,r3);
  m4=spline(rpr,mt4,r3);
  
  figure;
  loglog(r3,m1,'b',r3,m2,'r',r3,m3,'g',r3,m4,'m');
  saveas(gcf,sprintf('%s/%s_mprof_test1.png',result_dir,cluster));
  
  diff1=(m2-m1)./m1;
  diff2=(m3-m1)./m1;
  diff3=(m4-m1)./m1;
    
  figure;
  semilogx(r3,diff1,'r',r3,diff2,'g',r3,diff3,'m');
  saveas(gcf,sprintf('%s/%s_mprof_test2.png',result_dir,cluster));
  
  figure;
  loglog(r3,abs(diff1),'r',r3,abs(diff2),'g',r3,abs(diff3),'m');
  saveas(gcf,sprintf('%s/%s_mprof_test3.png',result_dir,cluster));
  
% figure;
% mm=mg+mst+mdm;
% loglog(rp,mt,'b',rp,mm.*1.1,'r');
% 
%rott=rot(rot>0);
%rpp=rp(rot>0);


%tpp=tp(tp>0);
%figure;
%loglog(rpp,rott)
%figure;
%loglog(rp,tp)

