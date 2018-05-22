function clustrules = init_clusterules()
clustrules.Compact          = 3; %  some cells with large dendrites could have values >2. I turned off filtering process in 
clustrules.diameter         = 10; % expected diameter of cells (used for 0.25 * pi/4*diam^2 < npixels < 30*pi/4*diam^2)

clustrules.npix_fraclow             = getOr(clustrules, {'npix_fraclow'}, (8/clustrules.diameter)^2);
clustrules.npix_frachigh            = getOr(clustrules, {'npix_frachigh'}, 50); %
