function cmap = colormap_symmetric()
% Create and activate blue-white-red colormap

n = 256;
half = n/2;

neg = [linspace(0,1,half)' linspace(0,1,half)' ones(half,1)];   % blue -> white
pos = [ones(half,1) linspace(1,0,half)' linspace(1,0,half)'];   % white -> red
cmap = [neg; pos];

colormap(cmap);
clim([-1,1] * max(abs(clim)));
