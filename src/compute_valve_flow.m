function ndot = compute_valve_flow(c_1, c_2, A ,C_d, T, M)
    % Compute molar flow through valve.
    % inputs:  c_1, c_2: gas concentration at left/right side in mol/m^3
    %          A: cross section in m^2
    %          C_d: drag coefficient
    %          T: temperature in K
    %          M: molar mass of components in kg/mol
    % output:  ndot: molar flow in mol/s

    % Universal gas constant in J/(mol*K)
    R = 8.314;

    % Pressure in Pa
    P_1 = sum(c_1) * R * T;
    P_2 = sum(c_2) * R * T;

    % Concentration of convective flow
    dP = P_1 - P_2;
    if dP > 0 
        c = c_1;
    else
        c = c_2;
    end
    
    % Density in kg/m^3 of convective flow
    rho = sum(c .* M);
    
    % Bernoulli: Mass flow in kg/s
    mdot = sign(dP)*C_d*A*sqrt(2*rho*abs(dP));
    
    % Fractions and molar mass
    y = c / sum(c);
    M_gas = sum(y .* M); % kg/mol

    % Total molar flow in mol/s
    ndot = y * mdot / M_gas;

end