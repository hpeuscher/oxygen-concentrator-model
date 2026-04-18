clearvars

params = devilbiss_default_params()

S = load('../datasets/20260410_122553_Devilbiss.mat'); % BMT2026 paper

history = S.history;

params.p_in = griddedInterpolant(history.t, history.p, 'linear', 'none');
params.n_comp = griddedInterpolant(history.t, history.n_des, 'linear', 'none');
params.f_patient = griddedInterpolant(history.t, history.f, 'linear', 'none');
params.psa_valve_state = S.psa_state;

tspan = [min(history.t), max(history.t)];

%% Run simulation
tic;
results = sim_psa(params, tspan, @get_bc_devilbiss);
toc;


%% Visualize results

plot_simulation_results_1d(results);
plot_simulation_results_2d(results, params);


figure; hold on; box on; grid on
plot(results.t,results.p_tank / 1e5,'LineWidth',2, 'DisplayName', 'simulation')
plot(history.t,history.p_outlet / 1e5,'LineWidth', 2, 'DisplayName', 'experiment')
xlabel('time in s'); ylabel('bar'); legend('show')
title('Pressure in tank')


figure; hold on; box on; grid on
plot(results.t,results.yO2_tank,'LineWidth',2, 'DisplayName', 'simulation');
plot(history.t,history.O2,'LineWidth',2, 'DisplayName', 'experiment');
xlabel('time in s'); ylabel('%'); legend('show')
title('O2 concentration at patient outlet')


figure; hold on; box on; grid on
plot(results.t, results.f_patient, 'LineWidth', 2, 'DisplayName', 'simulation');
plot(history.t,history.f, 'LineWidth', 2, 'DisplayName', 'experiment');
xlabel('time in s'); ylabel('L/min'); legend('show')
title('Patient Flow')


figure; hold on; box on; grid on;
plot(results.t, results.f_exhaust, 'LineWidth', 2, 'DisplayName', 'exhaust simulation');
plot(results.t, results.f_blow_1, 'LineWidth', 2, 'DisplayName', 'blowdown 1');
plot(results.t, results.f_blow_2, 'LineWidth', 2, 'DisplayName', 'blowdown 2');
plot(history.t, history.f_ex, 'LineWidth', 2, 'DisplayName', 'experiment');
xlabel('time in s'); ylabel('flow in L/min'); legend('show');
title("Flow at exhaust")


figure; hold on; box on; grid on;
plot(results.t, results.f_inlet, 'LineWidth', 2, 'DisplayName', 'simulation');
plot(history.t, history.f_in, 'LineWidth', 2, 'DisplayName', 'experiment');
xlabel('time in s'); ylabel('flow in L/min'); legend('show');
title("Flow from compressor")


figure; hold on; box on; grid on;
plot(history.t, history.n_des, 'LineWidth', 2);
xlabel('time in s'); ylabel('compressor speed in rpm'); legend('show');
title("Compressor speed")


% figure; hold on; box on; grid on;
% plot(results.t, results.power, 'LineWidth', 2, 'DisplayName', 'mechanical power');
% plot(hist.t, hist.power_el, 'LineWidth', 2, 'DisplayName', 'electrical power');
% xlabel('time in s'); ylabel('power in W'); legend('show');
% title("Power")

figure; hold on; box on; grid on;
stairs(history.t, history.psa_state, 'LineWidth', 2);
xlabel('time in s');
title("PSA State");


link_all_axes();

xlim(tspan([1,end]))
