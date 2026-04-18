function params = devilbiss_default_params()

N   = 30;                   % number of segments (dimensionless)

%% Natural constants and parameters
R = 8.314;                  % Universal gas constant in J/(mol*K)
V_mol = 22.7;               % molar volume in L/mol
P_atm = 101325;             % atmospheric pressure in Pa
D_L = 5e-5;                 % axial dispersion in m^2/s
mu = 1.85e-5;               % dynamic viscosity of air in Pa*s

T = 293;                    % temperature in K

%% Parameters for gas components (O2, N2)
y_atm = [0.21, 0.79];       % share of gas components in atmosphere
M = [0.032, 0.028];         % molar masses in kg/mol
k = [5, 5];                 % kinetic adsorption factor in 1/s
b = [0.2, 1.3];             % Langmuir coefficients in 1e-6/Pa
q_s = [2.5, 2.5];           % Adsorption capacities in mol/kg

%% Zeolite
eps = 0.3;                  % porosity (dimensionless)
dp = 1.6e-3;                % particle diameter in m
rho_p = 1250;               % particle density in kg/m^3

%% Geometry
L = 0.33;                   % cylinder length in m
D_col = 0.055;              % Cylinder diameter in m
A_cyl = pi * (D_col/2)^2;   % Cylinder cross section area in m^2
V_tank = 1e-3;              % Volume of tank in m^3
V_eq = 0.2e-3;              % Volume of equalization valve in m^3
V_exhaust = 1.5e-3;         % Volume of blowdown hose/exhaust in m^3

%% Compressor
comp_n_max = 7500;            % rpm
comp_n_min = 300;             % rpm
comp_p_max = 7;               % bar
comp_f_min = 75;              % L/min
comp_delta_f = 65;            % L/min

%% Valve cross section areas in mm^2
A_valve_out = 6;            % cylinder top -> tank
A_valve_bd = 24;            % cylinder bottom -> blowdown
A_valve_ex = 18;            % blowdown -> exhaust
A_valve_purge = 0.6;        % purge valve at top
A_valve_eq = 16;            % cylinder bottom <-> equalization volume
A_valve_patient = 4;        % tank -> patient outlet


%% Parameter structure
params = struct( ...
    "N", N, "L", L, "A_cyl", A_cyl, ...
    "V_tank", V_tank, "V_eq", V_eq, "V_exhaust", V_exhaust, ...
    "eps", eps, "dp", dp, "rho_p", rho_p, ...
    "comp_n_max", comp_n_max, "comp_n_min", comp_n_min, ...
    "comp_f_min", comp_f_min, "comp_delta_f", comp_delta_f,  "comp_p_max", comp_p_max, ...
    "R", R, "T", T, "V_mol", V_mol, "mu", mu, "D_L", D_L, ...
    "y_atm", y_atm, "P_atm", P_atm, "M", M, ...
    "b", b, "q_s", q_s, "k", k, ...
    "A_valve_out"    , A_valve_out    , ...
    "A_valve_bd"     , A_valve_bd     , ...
    "A_valve_ex"     , A_valve_ex     , ...
    "A_valve_eq"     , A_valve_eq     , ...
    "A_valve_purge"  , A_valve_purge  , ...
    "A_valve_patient", A_valve_patient ...
);
