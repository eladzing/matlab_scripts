%% explolre fb & xi dependence of stripping

%fg=0.1:0.1:0.3;%    0:0.01:0.25;
%beta=0.1:0.1:0.5; %       0.1:0.05:1;
fg=0.2;
beta=0.2;
BT=0:0.01:0.2;
fb=BT./(1-BT);
xi=1:0.5:10;

md=0.2;
mass=10.^10;
lambda=0.035;
len=length(fb)*length(xi);

ffb=reshape(repmat(fb,[length(xi) 1]),[len 1]);
xxi=repmat(xi,[1 length(fb)])';
mmass=mass.*ones(size(ffb));
ffg=fg.*ones(size(ffb));
bbeta=beta.*ones(size(ffb));

mmass=mmass./(1+ffg+ffb);  % set total mass as constant

if Flag
    tic
    res=rscale_mmw_array(mmass,'lambda',lambda,'md',md,...
        'fg',ffg,'beta',bbeta,'fb',ffb,'xi',xxi);
    toc
end
rd1=reshape(res.rd,[length(xi) length(fb)]);
si=mmass./(2*pi*res.rd.^2);
sig1=reshape(si,[length(xi) length(fb)]);

hf1=figure;
imagesc(BT,xi,log10(sig1))
set(gca,'Ydir','normal')
colorbar
xlabelmine('$B/T$')
ylabelmine('$\xi$')
titlemine('$\Sigma_s$')

%% calculate rps factor
pval=rps_factor_expdisk('etap',1,'sigma',mmass./(2*pi*res.rd.^2),'fd',ffg,...
    'Mc',1e15);
ppval=reshape(pval,[length(xi) length(fb)]);

hf2=figure;
imagesc(BT,xi,log10(ppval))
set(gca,'Ydir','normal')
colorbar
xlabelmine('$B/T$')
ylabelmine('$\xi$')
titlemine('$P_{rps}$')


r=0.01:0.01:100;
mStr=zeros(size(ffb));
for i=1:length(ffb)
    
    f=disk_force_reduced(r,'fg',ffg(i),'beta',bbeta(i),'fb',ffb(i),'xi',xxi(i));
    
    if max(f)>pval(i)
        rStr=interp1(f,r,pval(i));
        mStr(i)=exp_disk_mass(rStr,bbeta(i));
    else
        mStr(i)=0;
    end
    
end

mStr2d=reshape(mStr,[length(xi) length(fb)]);

hf3=figure;
imagesc(BT,xi,1-mStr2d)
set(gca,'Ydir','normal')
colorbar
xlabelmine('$B/T$')
ylabelmine('$\xi$')
titlemine('Mass stripped')
