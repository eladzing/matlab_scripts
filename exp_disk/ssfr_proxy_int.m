function ssfr=ssfr_proxy_int(ms,rd,fg,beta,fb) %,rmax)

%rmax=5.*rd;

knob=1;
A=(2.5+knob.*0.7).*1e-4;
alfa=1.4+knob.*0.15;
% sigma gas is in solar mass to pc^2
%siggas=mgas.*exp_disk_mass(rmax./rd,beta)./(pi.*(1e3.*rmax).^2); 

am1=alfa-1;

ssfr=(A*10^(-6*alfa)*(2*pi)^-am1/alfa^2).*(ms.^am1.*fg.^alfa.*beta.^(2*am1)./(rd.^(2*am1).*(1+fb))); 


% siggas=mgas./(pi.*(1e3.*rmax).^2); 
% @Article{Jaffe2015,
  Title                    = {{BUDHIES II: a phase-space view of H I gas stripping and star formation quenching in cluster galaxies}},
  Author                   = {{Jaff{\'e}}, Y.~L. and {Smith}, R. and {Candlish}, G.~N. and {Poggianti}, B.~M. and {Sheen}, Y.-K. and {Verheijen}, M.~A.~W. },
  Journal                  = {\mnras},
  Year                     = {2015},

  Month                    = apr,
  Pages                    = {1715-1728},
  Volume                   = {448},

  Abstract                 = {We investigate the effect of ram-pressure from the intracluster medium on the stripping of H I gas in galaxies in a massive, relaxed, X-ray bright, galaxy cluster at z = 0.2 from the Blind Ultra Deep H I Environmental Survey (BUDHIES). We use cosmological simulations, and velocity versus position phase-space diagrams to infer the orbital histories of the cluster galaxies. In particular, we embed a simple analytical description of ram-pressure stripping in the simulations to identify the regions in phase-space where galaxies are more likely to have been sufficiently stripped of their H I gas to fall below the detection limit of our survey. We find a striking agreement between the model predictions and the observed location of H I-detected and non-detected blue (late-type) galaxies in phase-space, strongly implying that ram-pressure plays a key role in the gas removal from galaxies, and that this can happen during their first infall into the cluster. However, we also find a significant number of gas-poor, red (early-type) galaxies in the infall region of the cluster that cannot easily be explained with our model of ram-pressure stripping alone. We discuss different possible additional mechanisms that could be at play, including the pre-processing of galaxies in their previous environment. Our results are strengthened by the distribution of galaxy colours (optical and UV) in phase-space, that suggests that after a (gas-rich) field galaxy falls into the cluster, it will lose its gas via ram-pressure stripping, and as it settles into the cluster, its star formation will decay until it is completely quenched. Finally, this work demonstrates the utility of phase-space diagrams to analyse the physical processes driving the evolution of cluster galaxies, in particular H I gas stripping.},
  Adsnote                  = {Provided by the SAO/NASA Astrophysics Data System},
  Adsurl                   = {http://adsabs.harvard.edu/abs/2015MNRAS.448.1715J},
  Archiveprefix            = {arXiv},
  Doi                      = {10.1093/mnras/stv100},
  Eprint                   = {1501.03819},
  Keywords                 = {galaxies: clusters: general, galaxies: clusters: invidivual: Abell 963, galaxies: clusters: intracluster medium, galaxies: evolution, galaxies: general},
  Owner                    = {eladzing},
  Timestamp                = {2015.02.23}
}

% sfr=(pi.*(rmax).^2).*A.*(siggas).^(1.4+0.15);
% 
% 
% %ssfr=sfr./ms.*exp_disk_mass(rmax./rd,1));
% ssfr=sfr./ms;

