list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%list2=[103 105 106 107];
result_dir='/home/titan3/eladzing/cold_flows/printout';

%pflag='print';
hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    
rmfac=0.02;rxfac=0.2;
bigbox=2;

for id=1:length(list1)
    halopath=sprintf(hpath,list1(id));
    clustername=sprintf('CL%d',list1(id))
    switch list1(id)
        case {103,105,106,107,5}
            smallbox=2
        otherwise
            smallbox=1
    end
    [full_ff full_ro full_ts rp m200 r200 t200]=catspheres(halopath,smallbox,bigbox);   
    ind=find( rp>rmfac.*r200 & rp<rxfac.*r200);
    innerr=ind(1)-1;
    outerr=ind(length(ind))+1;
    Mnorm=read_MGAS_Profile(halopath, r200)./1e13;
    
    yy=full_ff(innerr:outerr,:,:)./Mnorm;
    clear full_ts full_ro full_ff;
    y=yy(:); clear yy;
    n=hist(y,1000);
    dx=(max(y)-min(y))/1000;
    x=min(y):dx:max(y)-dx;
    figure;bar(x,log10(n));grid;
    title(sprintf( '%s Flux Histogram (0.02<r/R_{vir}<0.2)',clustername));
    xlabel('Mdot/Mg_{13}(R_{vir})');ylabel('log N');
    saveas(gcf,sprintf('%s/%s_fluxhist_core.png',result_dir,clustername));
    clear rp m200 r200 t200 x y n dx
    %clf;
end

