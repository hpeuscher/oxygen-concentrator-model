function plot_simulation_results_2d(results,params)

t = results.t;
N = size(results.cO2_1,2);
L = params.L;
z = linspace(0, L, N)';

%% ===================== GAS PHASE — CYLINDER 1 =====================
figure('Color', 'w', 'Position', [100 100 1000 600])

subplot(2,2,1)
surf(t, z, results.cO2_1', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
ylabel('z [m]')
title('Cyl. 1: c_{O_2} in mol m^{-3}')

subplot(2,2,2)
surf(t, z, results.cN2_1', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
title('Cyl. 1: c_{N_2} in mol m^{-3}')

subplot(2,2,3)
surf(t, z, results.yO2_1', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
xlabel('time t in s')
ylabel('z in m')
title('Cyl. 1: y_{O_2} in %')
caxis([0 100])

subplot(2,2,4)
surf(t, z, results.yN2_1', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
xlabel('time t in s')
title('Cyl. 1: y_{N_2} in %')
caxis([0 100])

%% ===================== GAS PHASE — CYLINDER 2 =====================
figure('Color', 'w', 'Position', [100 100 1000 600])

subplot(2,2,1)
surf(t, z, results.cO2_2', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
ylabel('z in m')
title('Cyl. 2: c_{O_2} in mol m^{-3}')

subplot(2,2,2)
surf(t, z, results.cN2_2', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
title('Cyl. 2: c_{N_2} in mol m^{-3}')

subplot(2,2,3)
surf(t, z, results.yO2_2', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
xlabel('time t in s')
ylabel('z [m]')
title('Cyl. 2: y_{O_2} in %')
caxis([0 100])

subplot(2,2,4)
surf(t, z, results.yN2_2', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
xlabel('time t in s')
title('Cyl. 2: y_{N_2} in %')
caxis([0 100])

%% ===================== GAS SPEED =====================
figure('Color', 'w', 'Position', [100 100 900 400])

subplot(1,2,1)
surf(t, z(1:end-1), results.u_1', 'EdgeAlpha', 0); view(2);
axis xy tight
colormap_symmetric(); colorbar();
xlabel('time t in s')
ylabel('z in m')
title('Cyl. 1: Gas Speed u(t,z) in m/s')

subplot(1,2,2)
surf(t, z(1:end-1), results.u_2', 'EdgeAlpha', 0); view(2);
axis xy tight
colormap_symmetric(); colorbar();
xlabel('time t in s')
title('Cyl. 2: Gas Speed u(t,z) in m/s')

%% ===================== ADSORPTION =====================
figure('Color', 'w', 'Position', [100 100 1000 600])

subplot(2,2,1)
surf(t, z, results.qO2_1', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
ylabel('z in m')
title('Cyl. 1: q_{O_2} in mol kg^{-1}')

subplot(2,2,2)
surf(t, z, results.qN2_1', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
title('Cyl. 1: q_{N_2} in mol kg^{-1}')

subplot(2,2,3)
surf(t, z, results.qO2_2', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
xlabel('Zeit t [s]')
ylabel('z in m')
title('Cyl. 2: q_{O_2} in mol kg^{-1}')

subplot(2,2,4)
surf(t, z, results.qN2_2', 'EdgeAlpha', 0); view(2);
axis xy tight
colorbar
xlabel('Zeit t [s]')
title('Cyl. 2: q_{N_2} in mol kg^{-1}')

%% ===================== PRESSURE =====================
figure('Color', 'w', 'Position', [100 100 900 400])

subplot(1,2,1)
surf(t, z, results.p_1' / 1e5, 'EdgeAlpha', 0); view(2);
axis xy tight
xlabel('time t in s')
ylabel('z in m')
colorbar
title('Cyl. 1: Pressure p(t,z) in bar')

subplot(1,2,2)
surf(t, z, results.p_2' / 1e5, 'EdgeAlpha', 0); view(2);
axis xy tight
xlabel('time t in s')
ylabel('z in m')
colorbar
title('Cyl. 2: Pressure p(t,z) in bar')

end
