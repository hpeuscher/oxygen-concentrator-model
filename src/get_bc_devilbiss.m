function BC = get_bc_devilbiss(t, c_1, c_2, c_tank, c_eqVol, c_bd, params)


C_d = 0.7;   % max valve flow coefficient
T = params.T;
M = params.M;

c_atm = params.y_atm * params.P_atm / (params.R * params.T);


A_valve_out     = 1e-6 * params.A_valve_out; % mm^2 -> m^2
A_valve_ex      = 1e-6 * params.A_valve_ex; % mm^2 -> m^2
A_valve_bd      = 1e-6 * params.A_valve_bd; % mm^2 -> m^2
A_valve_eq      = 1e-6 * params.A_valve_eq; % mm^2 -> m^2
A_valve_purge   = 1e-6 * params.A_valve_purge; % mm^2 -> m^2
A_valve_patient = 1e-6 * params.A_valve_patient; % mm^2 -> m^2


f_patient = params.f_patient;

[state, t_state] = params.psa_valve_state(t);
BC.tank = bc_tank_patient(c_tank);
BC.inlet_eq = bc_inlet_eq(c_eqVol);
BC.exhaust = bc_exhaust(c_bd);

BC.bottom1_ex = bc_cyl_bottom_blow(c_1(1,:), c_bd, state);
BC.bottom2_ex = bc_cyl_bottom_blow(c_2(1,:), c_bd, mod(state + 4, 8));
BC.bottom1_eq = bc_cyl_bottom_eq(c_1(1,:), c_eqVol, state);
BC.bottom2_eq = bc_cyl_bottom_eq(c_2(1,:), c_eqVol, mod(state + 4, 8));

BC.top1_tank = bc_cyl_top_tank(c_1(end,:), c_tank);
BC.top2_tank = bc_cyl_top_tank(c_2(end,:), c_tank);
BC.top1_top2 = bc_purge(c_1(end,:), c_2(end,:));


    function ndot = bc_inlet_eq(c_eqVol)
        % Flow from compressor to equalization volume
        nMax = params.comp_n_max;
        nMin = params.comp_n_min;
        fMin = params.comp_f_min;
        deltaf = params.comp_delta_f;
        deltaPMax = params.comp_p_max * 1e5; % bar -> Pa
        n = params.n_comp(t);
        P_eq = sum(c_eqVol) * (params.R * params.T);
        deltaP = P_eq - params.P_atm;
        Vdot_in = (n - nMin) / nMax * (fMin + (1 - deltaP / deltaPMax) * deltaf);

        if Vdot_in < 0
            Vdot_in = 0;
        end

        ndot = Vdot_in * params.y_atm / 60 / params.V_mol;
    end

    function ndot = bc_cyl_bottom_blow(c_cyl_bottom, c_bd, state)
        % Blowdown flow from cylinder bottom to exhaust
        switch state
            case 0 % adsorption
                C_d_blow = 0;
            case 1 % on the way to equalization
                C_d_blow = 0;
            case 2 % equalization and both cylinders pressurized
                C_d_blow = 0;
            case 3 % moving towards blowdown
                C_d_blow = C_d * sigmoid(t_state - 0.3, 0.2);
            case 4 % blowdown / purge
                C_d_blow = C_d;
            case 5 % on the way to equalization
                C_d_blow = C_d * (1 - sigmoid(t_state - 0.1, 0.2));
            case 6 % both cylinders pressurized
                C_d_blow = 0;
            case 7 % on the way to adsorption
                C_d_blow = 0;
        end

        ndot = compute_valve_flow(c_cyl_bottom, c_bd, A_valve_bd, C_d_blow, T, M);
    end

    function ndot = bc_cyl_bottom_eq(c_cyl_bottom, c_eqVol, state)
        % Equalization flow at cylinder bottom
        switch state
            case 0 % adsorption
                C_d_eq = C_d;
            case 1 % on the way to equalization
                C_d_eq = C_d;
            case 2 % equalization and both cylinders pressurized
                C_d_eq = C_d;
            case 3 % moving towards blowdown
                C_d_eq = C_d * (1 - sigmoid(t_state, 0.2));
            case 4 % blowdown / purge
                C_d_eq = 0;
            case 5 % on the way to equalization
                C_d_eq = C_d * sigmoid(t_state - 0.4, 0.1);
            case 6 % both cylinders pressurized
                C_d_eq = C_d;
            case 7 % on the way to adsorption
                C_d_eq = C_d;
        end

        ndot = compute_valve_flow(c_cyl_bottom, c_eqVol, A_valve_eq, C_d_eq, T, M);
    end

    function ndot = bc_cyl_top_tank(c_cyl_top, c_tank)
        % Flow from cylinder top to tank
        ndot = compute_check_valve_flow(c_cyl_top, c_tank, A_valve_out, C_d, T, M);
    end

    function ndot = bc_exhaust(c_bd)
        % Flow from blowdown hose to exhaust
        ndot = compute_check_valve_flow(c_bd, c_atm, A_valve_ex, C_d, T, M);
    end

    function ndot = bc_purge(c_top_1, c_top_2)
        % Purge flow between cylinder top
        ndot = compute_valve_flow(c_top_1, c_top_2, A_valve_purge, C_d, T, M);
    end

    function ndot = bc_tank_patient(c_tank)
        % Flow from tank to patient
        ndot = compute_check_valve_flow(c_tank, c_atm, A_valve_patient, C_d, T, M);
        if sum(ndot) > 0
            % max patient flow (limited by feedback controller)
            Vdot_max = f_patient(t); % target output flow in L/min
            Vdot = sum(ndot) * params.V_mol * 60; % in L/min

            % limiter
            ndot = ndot * min(1, Vdot_max / Vdot);
        end
    end
end

