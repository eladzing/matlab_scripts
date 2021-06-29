function re=find_re(fb,xi)

eta=0:0.01:10;
for i=1:length(fb)
ff=(1-exp(-eta).*(1+eta))+fb(i).*(xi(i).*eta).^2./(1+xi(i).*eta).^2;

re(i)=interp1(ff,eta,0.5);
end