function link_all_axes()
%% Link all x-axes

drawnow

figs = findall(0, 'Type', 'figure');
ax = findall(figs, 'Type', 'axes');
linkaxes(ax, 'x');
