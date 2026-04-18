function results = sim_psa(params, tspan, computeBCs)
% Simulate PSA process
% Inputs: params: parameter structure
%         tspan: simulation time span
%         computeBCs: function handle to compute boundary conditions


% Number of cells per cylinder
N = params.N;
% Segment length in m
dz  = params.L / N;

% Initial conditions...
P0 = params.P_atm;
c_atm0 = params.y_atm * P0 / (params.R * params.T);

% Initial concentrations in cylinder in mol/m^3
c_0 = ones(N,1)*c_atm0;

% Initial cylinder occupancy in mol/kg
q_0 = langmuir_equilibrium(c_0, params.b * 1e-6, kron(ones(N,1), params.q_s), params.R, params.T);

% Full initial state
x0 = [
    c_0(:)          % gas concentration in cylinder 1 in mol/m^3
    q_0(:)          % occupancy in cylinder 1 in mol/kg
    c_0(:)          % gas concentration in cylinder 2 in mol/m^3
    q_0(:)          % occupancy in cylinder 2 in mol/kg
    c_atm0(:)       % gas concentration in tank in mol/m^3
    c_atm0(:)       % gas concentration in equalization valve in mol/m^3
    c_atm0(:)       % gas concentration in exhaust volumne in mol/m^3
    ];

% Solve ODE system
[t,sol] = ode15s(@system_rhs, tspan, x0, odeset('AbsTol', 1e-3 * x0));

% Postprocess and return results
results = postprocess(t, sol);


    function dx = system_rhs(t,x)
        %% Returns right hand side of ODE system

        N = params.N;   % number of cells in each cylinder
        Ncomp = length(params.M);  % number of gas components

        % extract states of cylinder 1
        c_1 = reshape(x(1:Ncomp*N), N, Ncomp);
        q_1 = reshape(x(Ncomp*N+1:2*Ncomp*N), N, Ncomp);
        x(1:2*Ncomp*N) = [];

        % extract states of cylinder 2
        c_2 = reshape(x(1:Ncomp*N), N, Ncomp);
        q_2 = reshape(x(Ncomp*N+1:2*Ncomp*N), N, Ncomp);
        x(1:2*Ncomp*N) = [];

        % extract states of other volumes
        c_tank = x(1:2).';
        c_eq   = x(3:4).';
        c_exhaust = x(5:6).';

        % boundary conditions
        BCs = computeBCs(t, c_1, c_2, c_tank, c_eq, c_exhaust, params);

        % net cylinder flows
        ndot_bottom_1 = -BCs.bottom1_eq - BCs.bottom1_ex;
        ndot_bottom_2 = -BCs.bottom2_eq - BCs.bottom2_ex;
        ndot_top_1 = BCs.top1_tank + BCs.top1_top2;
        ndot_top_2 = BCs.top2_tank - BCs.top1_top2;

        % cylinder 1
        [dc_1, dq_1] = cylinder_rhs(c_1, q_1, ndot_bottom_1, ndot_top_1);

        % cylinder 2
        [dc_2, dq_2] = cylinder_rhs(c_2, q_2, ndot_bottom_2, ndot_top_2);

        % tank
        ndot_into_tank = BCs.top1_tank + BCs.top2_tank - BCs.tank;
        dc_tank = rhs_volume(ndot_into_tank, params.V_tank);

        % equalization volume
        ndot_into_eq = BCs.inlet_eq + BCs.bottom1_eq + BCs.bottom2_eq;
        dc_eq = rhs_volume(ndot_into_eq, params.V_eq);

        % blowdown hose / exhaust volume
        ndot_into_exhaust = BCs.bottom1_ex + BCs.bottom2_ex - BCs.exhaust;
        dc_exhaust = rhs_volume(ndot_into_exhaust, params.V_exhaust);

        % vector of state derivatives
        dx = [
            dc_1(:)
            dq_1(:)

            dc_2(:)
            dq_2(:)

            dc_tank(:)
            dc_eq(:)
            dc_exhaust(:)
            ];
    end

    function dc_dt = rhs_volume(ndot, V)
        %% Mass balance in control volume
        % inputs:  ndot: material flows in mol/s
        %          V: volume in m^3
        % output: concentration change in mol/(m^3*s)
        dc_dt = sum(ndot, 1) / V;
    end

    function [dc,dq] = cylinder_rhs(c, q, ndot_bottom, ndot_top)
        %% Mass balance in cylinder
        % inputs: c: concentration in cylinder cells in mol/m^3
        %         q: adsorbed material in cylinder cells in mol/kg
        %         ndot_bottom: inflow at bottom in mol/s
        %         ndot_top: outflot at top in mol/s
        % outputs: dc: change in concentration per cell in mol/(m^3*s)
        %          dq: change in adsorbed material in mol/(kg*s)

        % ------ adsorption ------
        % equilibrium of adsorption
        q_eq = langmuir_equilibrium(c, params.b * 1e-6, params.q_s, params.R, params.T);
        % LDF kinetics
        dq = (q_eq - q) * diag(params.k);

        % ------ convection ------
        c_left  = c(1:end-1,:);
        c_right = c(2:end,:);
        u_face = face_velocity(c);
        c_ip = (u_face >= 0) .* c_left + (u_face < 0) .* c_right;
        F_ip = u_face .* c_ip;
        conv_bottom = ndot_bottom / params.A_cyl;
        conv_top = ndot_top / params.A_cyl;
        % convective flow in m/s * mol/m^3
        convection = ([conv_bottom; F_ip] - [F_ip; conv_top]) / dz;

        % ------ dispersion ------
        % mirrored boundary conditions at the edges (no diffusion here)
        c_left  = [c(1,:); c(1:end-1,:)];
        c_right = [c(2:end,:); c(end,:)];
        % Fick's law
        diffusion = params.D_L * (c_left - 2*c + c_right) / dz^2;

        % ------ concentration change based on mass balance ------
        dc = (convection + diffusion - (1 - params.eps) * params.rho_p * dq) / params.eps;

    end

    function u_face = face_velocity(c)
        %% Compute face velocities using Ergun equation
        %   -(P(i+1) - P(i)) / dz = A*u + B*|u|*u

        % ----- cell-center pressure -----
        ctot = sum(c, 2);                     % mol/m^3
        P = ctot * params.R * params.T;       % Pa

        % ----- pressure drop -----
        dPdz_face = (P(2:end) - P(1:end-1)) / dz;   % (N-1)x1
        % absolute pressure gradient with mild regularization/smoothing
        dPdz_reg = hypot(dPdz_face, 1e-5);

        % ----- gas density at cell centers -----
        rho_g = c * params.M(:);              % kg/m^3

        % ----- face values for rho_g -----
        rho_face = 0.5 * (rho_g(1:end-1) + rho_g(2:end));   % (N-1)x1

        % ----- Ergun coefficients -----
        Ergun_coef = [150, 1.75];   % Ergun coefficients

        A = Ergun_coef(1) * params.mu * (1-params.eps)^2 ...
            / (params.dp^2 * params.eps^3);   % scalar

        B_face = Ergun_coef(2) * rho_face * (1-params.eps) ...
            / (params.dp * params.eps^3);     % (N-1)x1

        % ----- interior faces from Ergun -----
        % Solve:
        %   A*u + B*|u|*u = -(P(i+1)-P(i))/dz.
        u_face = sign(dPdz_face) .* (A - sqrt(A^2 + 4*B_face.*dPdz_reg)) ./ (2*B_face);

    end

    function q_eq = langmuir_equilibrium(c, b, q_s, R, T)
        %% Compute Langmuir equlibrium
        % inputs:  c: concentration of components in mol/m^3
        %          b: Langmuir coefficient of components in 1/Pa
        %          q_s: adsorption capacity of components in mol/kg
        % output:  q_eq: adsorbed mass in equilibrium in mol/kg

        % partial pressure (N x Ncomp) in Pa
        p = c * R * T;
        % equilibrial occupancy
        q_eq = q_s .* (p .* b) ./ (1 + p * b');
    end



    function res = postprocess(t, sol)
        %% Postprocessing
        res = table();
        res.t = t(:);
        N = params.N;
        R = params.R;
        T = params.T;

        % state variables sieve bed 1
        res.cO2_1 = sol(:, 1:N);
        res.cN2_1 = sol(:, 1*N+1:2*N);
        res.qO2_1 = sol(:, 2*N+1:3*N);
        res.qN2_1 = sol(:, 3*N+1:4*N);

        % state variables sieve bed 2
        res.cO2_2 = sol(:, 4*N+1:5*N);
        res.cN2_2 = sol(:, 5*N+1:6*N);
        res.qO2_2 = sol(:, 6*N+1:7*N);
        res.qN2_2 = sol(:, 7*N+1:8*N);

        % other state variables
        res.cO2_tank = sol(:,8*N+1);
        res.cN2_tank = sol(:,8*N+2);
        res.cO2_eq = sol(:,8*N+3);
        res.cN2_eq = sol(:,8*N+4);
        res.p_eq   = (res.cO2_eq + res.cN2_eq) * R*T;
        res.cO2_exhaust = sol(:,8*N+5);
        res.cN2_exhaust = sol(:,8*N+6);
        res.p_exhaust   = (res.cO2_exhaust + res.cN2_exhaust) * R*T;

        % derived quantities sieve bed 1
        res.cTot_1 = res.cO2_1 + res.cN2_1;
        res.yO2_1  = 100 * res.cO2_1 ./ res.cTot_1;
        res.yN2_1  = 100 * res.cN2_1 ./ res.cTot_1;
        res.p_1    = res.cTot_1 * R * T;

        % derived quantities sieve bed 2
        res.cTot_2 = res.cO2_2 + res.cN2_2;
        res.yO2_2  = 100 * res.cO2_2 ./ res.cTot_2;
        res.yN2_2  = 100 * res.cN2_2 ./ res.cTot_2;
        res.p_2    = res.cTot_2 * R * T;

        % --- boundary flow conditions  ---

        BCs = @(i) computeBCs(t(i), ...
            [res.cO2_1(i,:)', res.cN2_1(i,:)'], ...
            [res.cO2_2(i,:)', res.cN2_2(i,:)'], ...
            [res.cO2_tank(i), res.cN2_tank(i)], ...
            [res.cO2_eq(i), res.cN2_eq(i)], ...
            [res.cO2_exhaust(i), res.cN2_exhaust(i)], ...
            params);

        % bottom-side
        res.ndot_blow_1 = cell2mat(arrayfun(@(i) BCs(i).bottom1_ex(), (1:length(t))', 'UniformOutput', false));
        res.ndot_blow_2 = cell2mat(arrayfun(@(i) BCs(i).bottom2_ex(), (1:length(t))', 'UniformOutput', false));
        res.ndot_eq_1 = cell2mat(arrayfun(@(i) BCs(i).bottom1_eq(), (1:length(t))', 'UniformOutput', false));
        res.ndot_eq_2 = cell2mat(arrayfun(@(i) BCs(i).bottom2_eq(), (1:length(t))', 'UniformOutput', false));

        res.ndot_bottom_1 = -res.ndot_eq_1 - res.ndot_blow_1;
        res.ndot_bottom_2 = -res.ndot_eq_2 - res.ndot_blow_2;


        % equalization
        ndot_in_eq = arrayfun(@(i) BCs(i).inlet_eq(), (1:length(t))', 'UniformOutput', false);

        res.ndot_in_eq = cell2mat(ndot_in_eq);

        % top-side
        ndot_top1_tank = arrayfun(@(i) BCs(i).top1_tank(), (1:length(t))', 'UniformOutput', false);
        ndot_top2_tank = arrayfun(@(i) BCs(i).top2_tank(), (1:length(t))', 'UniformOutput', false);
        ndot_top1_top2 = arrayfun(@(i) BCs(i).top1_top2(), (1:length(t))', 'UniformOutput', false);
        res.ndot_top1_tank = cell2mat(ndot_top1_tank);
        res.ndot_top2_tank = cell2mat(ndot_top2_tank);
        res.ndot_top1_top2 = cell2mat(ndot_top1_top2);     
        res.ndot_top_1 = res.ndot_top1_tank + res.ndot_top1_top2;
        res.ndot_top_2 = res.ndot_top2_tank - res.ndot_top1_top2;

        % outflow to patient
        ndot_patient = arrayfun(@(i) BCs(i).tank(), (1:length(t))', 'UniformOutput', false);
        res.ndot_patient = cell2mat(ndot_patient);

        % outflow to exhaust
        ndot_exhaust = arrayfun(@(i) BCs(i).exhaust(), (1:length(t))', 'UniformOutput', false);
        res.ndot_exhaust = cell2mat(ndot_exhaust);

        % inputs
        res.p_in = params.p_in(t(:));


        % --- other derived quantities ---
        res.p_tank   = (res.cO2_tank + res.cN2_tank) * R * T;
        res.yO2_tank = 100 * res.cO2_tank ./ (res.cO2_tank + res.cN2_tank);
        res.yN2_tank = 100 * res.cN2_tank ./ (res.cO2_tank + res.cN2_tank);

        % flow speed
        res.u_1 = cell2mat(arrayfun(@(i)face_velocity([res.cO2_1(i,:)', res.cN2_1(i,:)'])', (1:length(t))', 'UniformOutput', false));
        res.u_2 = cell2mat(arrayfun(@(i)face_velocity([res.cO2_2(i,:)', res.cN2_2(i,:)'])', (1:length(t))', 'UniformOutput', false));

        % flow at inlet in L/min
        res.f_inlet = sum(res.ndot_in_eq, 2) * params.V_mol * 60;

        % flow at exhaust in L/min
        res.f_blow_1 = sum(res.ndot_blow_1, 2) * params.V_mol * 60;
        res.f_blow_2 = sum(res.ndot_blow_2, 2) * params.V_mol * 60;
        res.f_exhaust = sum(res.ndot_exhaust, 2) * params.V_mol * 60;

        % flow to patient in L/min
        res.f_patient = sum(res.ndot_patient, 2) * params.V_mol * 60;

        % mechanical power in W
        res.power = res.f_inlet / 1000 / 60  .* res.p_in;


        % total amount of substance in mol
        dV = params.A_cyl * params.L / params.N;
        rho = (1 - params.eps) * params.rho_p;
        res.nO2_1 = (sum(res.cO2_1, 2) * params.eps + sum(res.qO2_1, 2) * rho) * dV;
        res.nN2_1 = (sum(res.cN2_1, 2) * params.eps + sum(res.qN2_1, 2) * rho) * dV;
        res.nO2_2 = (sum(res.cO2_2, 2) * params.eps + sum(res.qO2_2, 2) * rho) * dV;
        res.nN2_2 = (sum(res.cN2_2, 2) * params.eps + sum(res.qN2_2, 2) * rho) * dV;

        % consistency check
        assert(max(abs(res.yO2_1(:) + res.yN2_1(:) - 100)) < 1e-10);
        assert(max(abs(res.yO2_2(:) + res.yN2_2(:) - 100)) < 1e-10);

    end


end
