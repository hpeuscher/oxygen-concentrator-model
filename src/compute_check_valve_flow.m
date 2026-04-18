function ndot = compute_check_valve_flow(c_1, c_2, A ,C_d, T, M)
% Compute molar flow through check valve from 1st to 2nd terminal.
% inputs:  c_1, c_2: gas concentration at terminals in mol/m^3
%          A: cross section in m^2
%          C_d: drag coefficient
%          T: temperature in K
%          M: molar mass of components in kg/mol
% output:  ndot: molar flow from 1st to 2nd terminal in mol/s

ndot = compute_valve_flow(c_1, c_2, A ,C_d, T, M);
if sum(ndot) < 0
    % check valve allows only outflow
    ndot = 0 * ndot;
end

end