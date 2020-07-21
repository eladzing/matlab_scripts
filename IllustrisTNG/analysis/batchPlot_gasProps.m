
printFlag=true;

entFlag=true;
tempFlag=false;
ssfrFlag=true;
bhratFlag=true;

for k=1:3
  
    switch k 
        case(1)
            gasField='Gal'
        case(2)
            gasField='CGM'
        case(3)
            gasField='Out'
    end
        
        
    
    plot_gasProperties_new
  
  readFlag=false;
  
end