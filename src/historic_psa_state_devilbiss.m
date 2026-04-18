function [phase, t_active] = historic_psa_state_devilbiss(state_t, state,t)
% Extracts the current PSA state from measurements
% inputs:  state_t:  timestamps of time series
%          state:    state time series (between 0 and 7)
%          t:        current timestamp
% outputs: phase:    current psa phase
%          t_active: time since start of current phase

idx = find(state_t >= t, 1, 'first');
if isempty(idx) || idx <= 1
    phase = 0;
    t_active = 0;
    return;
end
phase = state(idx - 1);
t_active = t - state_t(idx - 1);
end
