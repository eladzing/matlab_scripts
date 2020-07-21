list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

load mat_files/rho_profs2_a1.mat

srp=[];
srot=[];
srodm=[];

for ind=1:1
    
    %clname=rho_profs{ind,1};
    rp=rho_profs{ind,2};
    rot=rho_profs{ind,3};
    rodm=rho_profs{ind,4};
    virs=rho_profs{ind,5};
    rv=virs(1);
    mv=virs(2);
    
    norm=mv./(4.*pi.*rv.^3);
    srp(end+1,:)=log10(rp);
    srot(end+1,:)=log10(rot/norm);
    srodm(end+1,:)=log10(rodm/norm)-1;
end
    
 av_rot=mean(srot,1);
 std_rot=std(srot,0,1); 
 
 av_rod=mean(srodm,1);
 std_rod=std(srodm,0,1); 


