function [results, params] = sim_devilbiss_theoretical(params, tspan)

if nargin<1
    params = devilbiss_default_params();
end
if nargin<2
    tspan = [0, 120];
end

%% Inputs

if ~isfield(params, 'p_in')
    params.p_in = @(t) nan*t;
end
if ~isfield(params, 'n_comp')
    params.n_comp = @(t) 0*t + 7500;
end
if ~isfield(params, 'f_patient')
    params.f_patient = @(t) 0*t + 10;
end
if ~isfield(params, 'psa_valve_state')
    params.psa_valve_state = @get_psa_valve_state;
end



%% Run simulation
results = sim_psa(params, tspan, @get_bc_devilbiss);


%% Plot

plot_simulation_results_1d(results);
plot_simulation_results_2d(results, params);

figure; hold on; box on; grid on;
stairs(results.t, arrayfun(@(t)get_psa_valve_state(t), results.t), 'LineWidth', 2);
xlabel('time in s');
title("PSA State");

link_all_axes();
xlim(tspan([1,end]))


%% Determine current PSA phase
function [phase, t_state] = get_psa_valve_state(t)

t_ads = 2.48;
t_move = 0.5;
t_eq  = 0.1;

t_cycle = mod(t, 2*t_ads + 4*t_move + 2*t_eq);

if t_cycle < t_ads
    phase = 0;
    t_state = t_cycle;

elseif t_cycle < t_ads + t_move
    phase = 1;
    t_state = t_cycle - t_ads;

elseif t_cycle < t_ads + t_move + t_eq
    phase = 2;
    t_state = t_cycle - t_ads - t_move;

elseif t_cycle < t_ads + 2*t_move + t_eq
    phase = 3;
    t_state = t_cycle - t_ads - t_move - t_eq;

elseif t_cycle < 2*t_ads + 2*t_move + t_eq
    phase = 4;
    t_state = t_cycle - t_ads - 2*t_move - t_eq;

elseif t_cycle < 2*t_ads + 3*t_move + t_eq
    phase = 5;
    t_state = t_cycle - 2*t_ads - 2*t_move - t_eq;

elseif t_cycle < 2*t_ads + 3*t_move + 2*t_eq
    phase = 6;
    t_state = t_cycle - 2*t_ads - 3*t_move - t_eq;

else
    phase = 7;
    t_state = t_cycle - 2*t_ads - 3*t_move - 2*t_eq;

end
end
end