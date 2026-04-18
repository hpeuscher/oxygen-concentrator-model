function plot_simulation_results_1d(results)

t = results.t;


%% ===================== FIGURE - PRESSURE IN CYLINDER 1 =====================
figure; hold on; box on; grid on

plot(t,results.p_in / 1e5, 'k','LineWidth',2, 'DisplayName', 'Inlet');
plot(t,results.p_1(:,1) / 1e5,'b','LineWidth',2, 'DisplayName', 'Cylinder 1 bottom')
plot(t,results.p_eq / 1e5, 'y', 'LineWidth', 2, 'DisplayName', 'Equalization')
plot(t,results.p_1(:,end) / 1e5, 'r', 'LineWidth', 2, 'DisplayName', 'Cylinder 1 top')
plot(t,results.p_tank / 1e5, 'g', 'LineWidth', 2,'DisplayName', 'Tank');
plot(t,results.p_exhaust / 1e5, 'DisplayName', 'Exhaust');
legend show;
title('Pressure in Cylinder 1')
ylabel('bar')
xlabel('Time in s')

%% ===================== FIGURE - PRESSURE IN CYLINDER 2 =====================
figure; hold on; box on; grid on
plot(t,results.p_in / 1e5, 'k', 'LineWidth', 2,'DisplayName', 'Inlet');
plot(t,results.p_2(:,1) / 1e5, 'b', 'LineWidth', 2, 'DisplayName', 'Cylinder 2 bottom')
plot(t,results.p_eq / 1e5, 'y', 'LineWidth', 2, 'DisplayName', 'Equalization')
plot(t,results.p_2(:,end) / 1e5, 'r', 'LineWidth', 2, 'DisplayName', 'Cylinder 2 top')
plot(t,results.p_tank / 1e5, 'g', 'LineWidth', 2, 'DisplayName', 'Tank');
plot(t,results.p_exhaust / 1e5, 'DisplayName', 'Exhaust');
legend show;
title('Pressure in Cylinder 2')
ylabel('bar')
xlabel('Time in s')

%% ===================== FIGURE - Tank =====================
figure;

subplot(3,1,1); hold on; box on; grid on
plot(t,results.p_tank / 1e5,'LineWidth', 2);
ylabel('bar')
title('Tank pressure')

subplot(3,1,2); hold on; box on; grid on
plot(t,results.yO2_tank, 'LineWidth', 2);
ylabel('%')
title('Fraction of Oxygen in Tank')

subplot(3,1,3); hold on; box on; grid on
plot(t,results.yN2_tank, 'LineWidth', 2);
ylabel('%')
xlabel('Time in s')
title('Fraction of Nitrogen in Tank')


%% ===================== FIGURE - EQ =====================
figure; hold on; box on; grid on
plot(t, sum(results.ndot_in_eq, 2), 'LineWidth', 2, 'DisplayName', 'inlet');
plot(t, sum(results.ndot_eq_1, 2), 'LineWidth', 2, 'DisplayName', 'cylinder 1');
plot(t, sum(results.ndot_eq_2, 2), 'LineWidth', 2, 'DisplayName', 'cylinder 2');
ylabel('mol/s')
title('Flows at equalization valve')
xlabel('Time in s')
legend show;

%% ===================== FIGURE - CYLINDER TOP =====================
figure; hold on; box on; grid on
plot(t, sum(results.ndot_top1_tank, 2), 'LineWidth', 2, 'DisplayName', 'top1->tank');
plot(t, sum(results.ndot_top2_tank, 2), 'LineWidth', 2, 'DisplayName', 'top2->tank');
plot(t, sum(results.ndot_top1_top2, 2), 'LineWidth', 2, 'DisplayName', 'top1->top2');
ylabel('mol/s')
title('Flows at cylinder top')
xlabel('Time in s')
legend show;


end
