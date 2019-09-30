function [ tout, recording_vars] = Integration(inputs,mode,tspan)
% Integrate necesary differential equations for rocket engine modeling 
% using Euler's method. This includes real combustion properties and 
% supercharging.

%   INPUTS:
%       - inputs: structure of motor characteristics (all units SI: m, s, 
%       kg, K, mol)
%           - CombustionData: string filename in "Combustion Data/" folder 
%           from which to source combustion data
%           - Injexit_area: total orifice area of injector
%           - Cdischarge: discharge coefficient of injector orifices
%           - MOVTime: time for injector flow to ramp up to 100%
%           - rocket_dry_mass: dry mass of rocket
%           - tankvol: volume of oxidizer tank
%           - l_vol: initial volume of liquid nitrous oxide in oxidizer 
%           tank
%           - tank_id: inner diameter of tank
%           - h_offset: height difference between bottom of tank and 
%           injector
%           - flowline_id: inner diameter of flow line
%           - Ttank: tank initial temperature
%           - Fueldensity: density of fuel
%           - grainlength: length of fuel grain
%           - graindiameter: outer diameter of grain
%           - portradius0: grain intial port radius
%           - chamberlength: combustion chamber total length
%           - n: grain ballistic coefficient (r_dot = a*G^n)
%             - a: grain ballistic coefficient (r_dot = a*G^n)
%           - nozzle_efficiency: nozzle exhaust energy efficiency 
%           (proportion)
%           - nozzle_correction_factor: proportion of ideal thrust 
%           (factor for divergence losses, etc.) actually achieved
%           - c_star_efficiency: combustion efficiency / proportion of c*
%           actually achieved
%           - Tdiameter: nozzle throat diameter
%           - E: expansion ratio
%           - Tamb: ambient temperature
%           - Pamb: ambient pressure
%           - SPress: supercharging regulator pressure
%           - M_sc: molecular mass
%           - SVol: volume of external pressurization tank (0 for 
%           supercharging / no external tank)
%           - c_v_S: specific heat at constant volume of pressurant gas
%             - c_p_S: specific heat at constant volume of pressurant gas
%           - P_S_i: initial pressurant gas storage pressure (must be 
%           present, but not used for supercharging / no external tank)
%           - S_CdA: effective flow area of pressurant gas (must be 
%           but not used for supercharging / no external tank)
%           - Charged: 1 for pressurant gas present, 0 for no pressurant 
%           gas present
%       - mode: structure of options defining mode of motor operation
%           - combustion_on: 1 for hot-fire, 0 for cold-flow
%           - flight_on: 1 for flight conditions, 0 for ground conditions
%       - tspan: vector of time values over which to record outputs
%   OUTPUS:
%       - tspan: output time vector
%       - F_thrust: thrust over tspan
%       - p_cc: combustion chamber pressure over tspan
%       - p_oxtank: tank pressure over tspan
%       - p_oxmanifold: oxidizer manifold pressure over tspan
%       - T_tank: tank temperature over tspan
%       - T_cc: combustion chamber as temperature over tspan
%       - area_core: fuel grain core area over tspan
%       - OF: oxidizer/fuel ratio over tspan
%       - m_dot_ox: oxidizer mass flow rate over tspan
%       - m_dot_ox_crit: critical two-phase oxidizer mass flow rate over tspan
%       - p_crit: two-phase critical downstream pressure over tspan
%       - M_e: exit Mach number over tspan
%       - p_exit: nozzle exit pressure over tspan
%       - p_shock: critical back pressure for normal shock formation over
%       tspan

%% State Vector Initialization
if mode.type == 'hybrid'
    [state_0, x0] = InitializeHybridState(inputs, mode);
elseif mode.type == 'liquid'
    [state_0, x0] = InitializeLiquidState(inputs, mode);
end

%% Solution
InitializeRecord(tspan, {'F_thrust', 'p_cc', 'p_oxtank', 'p_oxpresstank', ...
    'p_fueltank', 'p_fuelpresstank', 'p_oxmanifold', 'T_oxtank', 'T_cc', ...
    'area_core', 'gamma_ex', 'm_lox', 'm_gox', 'm_fuel', 'p_crit', ...
    'm_dot_ox_crit', 'M_e', 'p_exit', 'p_shock'});
options = odeset('Events', @(t,y) TerminalFunction(t,y,inputs, mode),...
    'MaxStep',0.01);
if strcmp(mode.type,'liquid')
    odefun = @(t,y) LiquidModel(t,y,inputs,mode);
elseif strcmp(mode.type,'hybrid')
    odefun = @(t,y) HybridModel(t,y,inputs,mode);
end
[tout,x,te,xe,ie] = ode15s(odefun,tspan,x0,options);

%% Outputs
recording_vars = OutputRecord(tout, inputs, mode);
recording_vars.m_press = 0;
if inputs.ox_pressurant.active
    recording_vars.m_press = recording_vars.m_press + ...
        state_0.m_oxtank_press + state_0.m_oxpresstank;
end
if strcmp(mode.type, 'liquid') && inputs.fuel_pressurant.active
    recording_vars.m_press = recording_vars.m_press + ...
        state_0.m_fueltank_press + state_0.m_fuelpresstank;
end

% Include event times in record
if any(ie == 2)
    recording_vars.t_fuel_liq = te(ie==2);
else
    recording_vars.t_fuel_liq = NaN;
end
if any(ie == 4)
    recording_vars.t_ox_liq = te(ie == 4);
else
    recording_vars.t_ox_liq = NaN;
end
end

function x_dot = HybridModel(time,x,inputs,mode)
% Model engine physics for a hybrid motor, provide state vector derivative

%% Initialization of State Properties
state = HybridStateVector.FromColumnVector(x, inputs, mode);

state_dot = HybridStateVector();

%% Calculate Injector Mass Flow Rate
[m_dot_lox, m_dot_gox, m_dot_oxtank_press, T_dot_drain, p_crit, m_dot_ox_crit] = ...
    N2OTankMDot(inputs, state, time);
[m_dot_vap, T_dot_vap] = N2OTankEquilibrium(inputs, state);

state_dot.m_lox = -m_dot_lox - m_dot_vap;
state_dot.m_gox = -m_dot_gox + m_dot_vap;
state_dot.m_oxtank_press = -m_dot_oxtank_press;
state_dot.T_oxtank = T_dot_drain + T_dot_vap; 

m_dot_ox = m_dot_lox + m_dot_gox;

%% Combustion Dynamics
if mode.combustion_on % for hot fire only 
    [m_dot_f, r_dot] = HybridMdotFuel(inputs, state, m_dot_ox);
    
    [F_thrust, OF, M_e, p_exit, p_shock, m_cc_dot, M_cc_dot, ...
    gamma_cc_dot, T_cc_dot] = CombustionChamber(inputs, state, m_dot_ox, m_dot_f);
else % no combustion (cold flow)
    m_dot_f = 0;
    r_dot = 0;
    
    F_thrust = 0;
    OF = 0;
    M_e = 0;
    p_exit = inputs.p_amb;
    p_shock = 0;
    m_cc_dot = 0;
    M_cc_dot = 0;
    gamma_cc_dot = 0;
    T_cc_dot = 0;
end
state_dot.m_fuel = -m_dot_f;
state_dot.m_cc = m_cc_dot;
state_dot.M_cc = M_cc_dot;
state_dot.gamma_cc = gamma_cc_dot;
state_dot.T_cc = T_cc_dot;

%% Pressurant Flow
if inputs.ox_pressurant.active
    [T_dot_press, m_dot_press] = PressurantFlow(state.p_oxpresstank, ...
        state.T_oxpresstank, state.p_oxtank, state.T_oxtank, ...
        state.oxtank_m_cv, inputs.ox_pressurant);
else
    T_dot_press = 0;
    m_dot_press = 0;
end
state_dot.T_oxtank = state_dot.T_oxtank + T_dot_press;
state_dot.m_oxtank_press = m_dot_press;
state_dot.m_oxpresstank = -m_dot_press;

%% Generate column vector
x_dot = state_dot.ColumnVector;
HybridRecord(time, state, state_dot, F_thrust, p_crit, ...
    m_dot_ox_crit, M_e, p_exit, p_shock);

% fprintf('time: %.3f \n', time)
% fprintf('%.3g\t', x_dot)
% fprintf('\n')
% fprintf('%% (Ptank - Pvap)/Pvap: %.1f\n', 100*(state.p_gox - state.N2O_properties.Pvap)/state.N2O_properties.Pvap)

end

function x_dot = LiquidModel(time,x,inputs,mode)
% Model engine physics for a liquid motor, provide state vector derivative

%% Constants
R_u = 8.3144621; % Universal gas constant [J/mol*K]
M_n2o = 0.044013; % Molecular mass of nitrous oxide [kg/mol]
R_n2o = R_u/M_n2o; %  Specific gas constant of nitrous oxide [J/kg*K]
a_n2o = 0.38828/M_n2o^2; % van der Waal's constant a for N2O [[Pa*(kg/m^3)^-2]]
b_n2o = 44.15/M_n2o*10^-6; % van der Waal's constant a for N2O [m^3/kg]

%% Initialization of State Properties
state = LiquidStateVector.FromColumnVector(x, inputs, mode);

state_dot = LiquidStateVector();

%% Calculate Injector Mass Flow Rate
[m_dot_lox, m_dot_gox, m_dot_oxtank_press, T_dot_drain_ox, p_crit, m_dot_ox_crit] = ...
    N2OTankMDot(inputs, state, time);

[m_dot_f, T_dot_drain_f] = FuelTankMDot(inputs, state, time);

[m_dot_vap, T_dot_vap] = N2OTankEquilibrium(inputs, state);

state_dot.m_lox = -m_dot_lox - m_dot_vap;
state_dot.m_gox = -m_dot_gox + m_dot_vap;
state_dot.m_oxtank_press = -m_dot_oxtank_press;
state_dot.T_oxtank = T_dot_drain_ox + T_dot_vap; 

state_dot.m_fuel = -m_dot_f;
state_dot.T_fueltank_press = T_dot_drain_f; 

m_dot_ox = m_dot_lox + m_dot_gox;

%% Combustion Dynamics
if mode.combustion_on % for hot fire only 
    [F_thrust, OF, M_e, p_exit, p_shock, m_cc_dot, M_cc_dot, ...
    gamma_cc_dot, T_cc_dot] = CombustionChamber(inputs, state, m_dot_ox, m_dot_f);
else % no combustion (cold flow)
    F_thrust = 0;
    OF = 0;
    M_e = 0;
    p_exit = inputs.p_amb;
    p_shock = 0;
    m_cc_dot = 0;
    M_cc_dot = 0;
    gamma_cc_dot = 0;
    T_cc_dot = 0;
end
state_dot.m_cc = m_cc_dot;
state_dot.M_cc = M_cc_dot;
state_dot.gamma_cc = gamma_cc_dot;
state_dot.T_cc = T_cc_dot;

%% Oxidizer Pressurant Flow
if inputs.ox_pressurant.active
    [T_dot_press, m_dot_press] = PressurantFlow(state.p_oxpresstank, ...
        state.T_oxpresstank, state.p_oxtank, state.T_oxtank, ...
        state.oxtank_m_cv, inputs.ox_pressurant);
else
    T_dot_press = 0;
    m_dot_press = 0;
end
state_dot.T_oxtank = state_dot.T_oxtank + T_dot_press;
state_dot.m_oxtank_press = m_dot_press;
state_dot.m_oxpresstank = -m_dot_press;

%% Fuel Pressurant Flow
if inputs.fuel_pressurant.active
    [T_dot_press, m_dot_press] = PressurantFlow(state.p_fuelpresstank, ...
        state.T_fuelpresstank, state.p_fueltank, state.T_fueltank_press, ...
        state.m_fueltank_press*inputs.fuel_pressurant.gas_properties.c_v, ...
        inputs.fuel_pressurant);
else
    T_dot_press = 0;
    m_dot_press = 0;
end
state_dot.T_fueltank_press = state_dot.T_fueltank_press + T_dot_press;
state_dot.m_fueltank_press = m_dot_press;
state_dot.m_fuelpresstank = -m_dot_press;

%% Generate column vector
x_dot = state_dot.ColumnVector;
LiquidRecord(time, state, state_dot, F_thrust, p_crit, ...
    m_dot_ox_crit, M_e, p_exit, p_shock);

% fprintf('time: %.3f \n', time)
% fprintf('%.3g\t', x_dot)
% fprintf('\n')
% fprintf('%% (Ptank - Pvap)/Pvap: %.1f\n', 100*(state.p_gox - state.N2O_properties.Pvap)/state.N2O_properties.Pvap)

end

function [value,isterminal,direction] = TerminalFunction(t,x,inputs, mode)
%TerminalFunction Sets the termination of integration. Integration will 
%stop when values go to zero. 

if strcmp(mode.type, 'liquid')
    state = LiquidStateVector.FromColumnVector(x, inputs, mode);
elseif strcmp(mode.type, 'hybrid')
    state = HybridStateVector.FromColumnVector(x, inputs, mode);
end

isterminal = [1, 1];
direction = [0, 0];

value = [state.m_lox + state.m_gox, state.m_fuel];  % mass of propellants
end

function [state_0, x0] = InitializeHybridState(inputs, mode)
%INITIALIZEHYBRIDSTATE Initializes the state vector for a hybrid system.
%   Uses the inputs to create an initial state vector for a hybrid

% Initialize State
state_0 = HybridStateVector(inputs, mode);

state_0 = InitializeOxtank(state_0, inputs);

state_0.m_fuel = pi*((inputs.fuel.grain_od/2)^2 - inputs.fuel.port_rad^2)*...
            inputs.fuel.grain_length*inputs.fuel.rho;
        
state_0 = InitializeCombustionChamber(state_0, inputs);

if state_0.V_cc <= 0
    error('Grain exceeds combustion chamber volume.')
end

x0 = state_0.ColumnVector;
end

function [state_0, x0] = InitializeLiquidState(inputs, mode)
%INITIALIZELIQUIDSTATE Initializes the state vector for a liquid system.
%   Uses the inputs to create an initial state vector for a liquid

% Initialize State
state_0 = LiquidStateVector(inputs, mode);

state_0 = InitializeOxtank(state_0, inputs);

state_0 = InitializeFueltank(state_0, inputs);
        
state_0 = InitializeCombustionChamber(state_0, inputs);

x0 = state_0.ColumnVector;
end

function [state] = InitializeOxtank(state, inputs)
% Initialize oxidizer tank based on inputs

state.T_oxtank = inputs.ox.T_tank;
state.N2O_properties = N2O_Properties(state.T_oxtank);
state.m_lox = state.N2O_properties.rho_l*inputs.ox.V_l;
state.m_gox = state.N2O_properties.rho_g*state.V_ox_ullage;
if inputs.ox_pressurant.active
    state.m_oxtank_press = (inputs.ox_pressurant.set_pressure - state.p_gox)*...
        state.V_ox_ullage/(inputs.ox_pressurant.gas_properties.R_specific*state.T_oxtank);
    state.m_oxpresstank = inputs.ox_pressurant.storage_initial_pressure*...
        inputs.ox_pressurant.tank_volume/(inputs.ox_pressurant.gas_properties.R_specific*inputs.T_amb);
else
    state.m_oxtank_press = 0;
    state.m_oxpresstank = 0;
end

end

function [state] = InitializeFueltank(state, inputs)
% Initialize oxidizer tank based on inputs

state.T_fueltank_press = inputs.T_amb;
state.m_fuel = inputs.fuel.V_l*inputs.fuel.rho;
if inputs.fuel_pressurant.active
    state.m_fueltank_press = inputs.fuel_pressurant.set_pressure*...
        state.V_fuel_ullage/(inputs.fuel_pressurant.gas_properties.R_specific*state.T_fueltank_press);
    state.m_fuelpresstank = inputs.fuel_pressurant.storage_initial_pressure*...
        inputs.fuel_pressurant.tank_volume/(inputs.fuel_pressurant.gas_properties.R_specific*inputs.T_amb);
else
    state.m_fueltank_press = 0;
    state.m_fuelpresstank = 0;
end

end

function [state] = InitializeCombustionChamber(state, inputs)
% Initialize combustion chamber based on inputs

% Constants
M_air = 0.02897; % mean molecular mass of air [kg/mol]
rho_air = 1.225; % kg/m^3
gamma_air = 1.4;

state.m_cc = state.V_cc*rho_air;
state.M_cc = M_air;
state.gamma_cc = gamma_air;
state.T_cc = inputs.T_amb;
end

function InitializeRecord(time, recording_vars)
%InitializeRecord Initialize the recording variables.

global recording_values dt n_rec

n_rec = 0;
recording_values = cell2struct(cell(length(recording_vars),1),recording_vars);

dt = mean(diff(time));
end

function HybridRecord(time, state, state_dot, F_thrust, p_crit, ...
    m_dot_ox_crit, M_e, p_exit, p_shock)
% Record variables not recorded by state

global recording_values dt n_rec time_rec

if floor(time/dt) > n_rec && (n_rec == 0 || time > time_rec(n_rec))
    time_rec(n_rec + 1) = time;
    recording_values.F_thrust(n_rec + 1) = F_thrust;
    recording_values.p_cc(n_rec + 1) = state.p_cc;
    recording_values.p_oxtank(n_rec + 1) = state.p_oxtank;
    recording_values.p_oxpresstank(n_rec + 1) = state.p_oxpresstank;
    recording_values.p_fueltank(n_rec + 1) = NaN;
    recording_values.p_fuelpresstank(n_rec + 1) = NaN;
    recording_values.p_oxmanifold(n_rec + 1) = state.p_oxmanifold;
    recording_values.T_oxtank(n_rec + 1) = state.T_oxtank;
    recording_values.T_cc(n_rec + 1) = state.T_cc;
    recording_values.m_fuel(n_rec + 1) = state.m_fuel;
    recording_values.area_core(n_rec + 1) = pi*state.rad_port^2;
    recording_values.gamma_ex(n_rec + 1) = state.gamma_cc;
    recording_values.m_lox(n_rec + 1) = state.m_lox;
    recording_values.m_gox(n_rec + 1) = state.m_gox;
    recording_values.p_crit(n_rec + 1) = p_crit;
    recording_values.m_dot_ox_crit(n_rec + 1) = m_dot_ox_crit;
    recording_values.M_e(n_rec + 1) = M_e;
    recording_values.p_exit(n_rec + 1) = p_exit;
    recording_values.p_shock(n_rec + 1) = p_shock;

    n_rec = n_rec + 1;
end
end

function LiquidRecord(time, state, state_dot, F_thrust, p_crit, ...
    m_dot_ox_crit, M_e, p_exit, p_shock)
% Record variables not recorded by state

global recording_values dt n_rec time_rec

if floor(time/dt) > n_rec && (n_rec == 0 || time > time_rec(n_rec))
    time_rec(n_rec + 1) = time;
    recording_values.F_thrust(n_rec + 1) = F_thrust;
    recording_values.p_cc(n_rec + 1) = state.p_cc;
    recording_values.p_oxtank(n_rec + 1) = state.p_oxtank;
    recording_values.p_oxpresstank(n_rec + 1) = state.p_oxpresstank;
    recording_values.p_fueltank(n_rec + 1) = state.p_fueltank;
    recording_values.p_fuelpresstank(n_rec + 1) = state.p_fuelpresstank;
    recording_values.p_oxmanifold(n_rec + 1) = state.p_oxmanifold;
    recording_values.T_oxtank(n_rec + 1) = state.T_oxtank;
    recording_values.T_cc(n_rec + 1) = state.T_cc;
    recording_values.m_fuel(n_rec + 1) = state.m_fuel;
    recording_values.area_core(n_rec + 1) = NaN;
    recording_values.gamma_ex(n_rec + 1) = state.gamma_cc;
    recording_values.m_lox(n_rec + 1) = state.m_lox;
    recording_values.m_gox(n_rec + 1) = state.m_gox;
    recording_values.p_crit(n_rec + 1) = p_crit;
    recording_values.m_dot_ox_crit(n_rec + 1) = m_dot_ox_crit;
    recording_values.M_e(n_rec + 1) = M_e;
    recording_values.p_exit(n_rec + 1) = p_exit;
    recording_values.p_shock(n_rec + 1) = p_shock;

    n_rec = n_rec + 1;
end
end

function [record] = OutputRecord(tspan, inputs, mode)
% Output values of recorded variables

global recording_values n_rec time_rec

record = recording_values;

time_rec(n_rec + 1:end) = [];

fields = fieldnames(recording_values);
for ii = 1:length(fields)
    record.(fields{ii})(n_rec + 1:end) = [];
    record.(fields{ii}) = interp1(time_rec, ...
        record.(fields{ii}), tspan, 'linear', 'extrap');
end

A_star = (pi/4*inputs.d_throat^2);
record.time = tspan;
record.m_ox = record.m_lox + record.m_gox;
record.impulse = trapz(tspan,record.F_thrust);
record.m_dot_ox = -[0; diff(record.m_ox)./diff(tspan)];
record.m_dot_fuel = -[0; diff(record.m_fuel)./diff(tspan)];
record.OF_i = record.m_dot_ox./record.m_dot_fuel;
record.OF = trapz(tspan,record.m_dot_ox)/trapz(tspan,record.m_dot_fuel);
record.m_dot_prop = record.m_dot_ox + record.m_dot_fuel;
record.c_star_i = record.p_cc*A_star./record.m_dot_prop;
record.c_star = trapz(tspan,record.p_cc)*A_star...
    ./trapz(tspan,record.m_dot_prop);
record.c_f_i = record.F_thrust./(A_star*record.p_cc);
record.c_f = record.impulse...
    ./(trapz(tspan,record.p_cc)*A_star);
record.Isp_i = record.F_thrust./record.m_dot_prop;
record.Isp = record.impulse./trapz(tspan,record.m_dot_prop);

record.ox_pressure_drop = (record.p_oxmanifold-record.p_cc)./record.p_oxtank;
if mode.type == 'liquid'
    record.fuel_pressure_drop = (record.p_fueltank-record.p_cc)./record.p_fueltank;
end

clear recording_values n_rec dt time_rec
end

function [m_dot_lox, m_dot_gox, m_dot_oxtank_press, T_dot_drain, p_crit, m_dot_ox_crit] = ...
    N2OTankMDot(inputs, state, time)
%N2OTankMDot Calculates flow rate out of N2O Tank
%   Interpolates between liquid and gas flow based on amount of liquid
%   oxidizer remaining. If no liquid oxidizer, uses gas flow. If plenty of
%   liquid oxidizer (i.e. mass greater than set tolerance), liquid flow is
%   used. Linear interpolation in between (in order to avoid hysteresis
%   with a small, steady liquid source in tank. 
%
%   Inputs:
%       - inputs: object representing motor inputs
%       - state: object representing state of system
%       - time: time
%   Outputs:
%       - m_dot_lox: mass flow rate out of tank of liquid oxidizer
%       - m_dot_gox: mass flow rate out of tank of gaseous oxidizer
%       - m_dot_oxtank_press: mass flow rate out of tank of pressurant gas
%       - T_dot_drain: oxidizer tank temperature change from draining
%       - p_crit: critical pressure below which draining flow is choked
%       - m_dot_ox_crit: critical mass flow rate below which draining flow 
%         is choked

% Options
global dm_lox_tol
dm_lox_tol = 1e-2; % kg

% Constants
M_n2o = 0.044013; % Molecular mass of nitrous oxide [kg/mol]
a_n2o = 0.38828/M_n2o^2; % van der Waal's constant a for N2O [[Pa*(kg/m^3)^-2]]

% If tank pressure is greater than combustion chamber pressure
if (state.p_oxmanifold > state.p_cc)&&(inputs.ox.Cd_injector*inputs.ox.injector_area*inputs.Throttle(time)>0)
    % Find flow rates for gas, flow rates for liquid, interpolate within
    % tolerance for smooth transition
    m_dot_lox = zeros(2,1);
    m_dot_gox = zeros(2,1);
    m_dot_oxtank_press = zeros(2,1);
    p_crit = zeros(2,1);
    m_dot_ox_crit = zeros(2,1);
    Q = zeros(2,1);
    
    %% Liquid flow
    [m_dot_ox, m_dot_ox_crit(1), p_crit(1) ]= LN2OMDot(inputs, ...
        inputs.ox.injector_area*inputs.Throttle(time), state.p_oxmanifold, state.T_oxtank, ...
        state.p_cc);

    m_dot_lox(1) = m_dot_ox;
    m_dot_gox(1) = 0;
    m_dot_oxtank_press(1) = 0;
    Q(1) = m_dot_ox/state.N2O_properties.rho_l;
    
    %% Gas Flow
    d_inj = sqrt(4/pi*inputs.ox.injector_area*inputs.Throttle(time));
    [ ~, ~, ~, ~, m_dot_ox, ~ ] = NozzleCalc( d_inj, d_inj, state.T_oxtank, ...
        state.p_oxmanifold, state.gamma_ox_ullage, M_n2o, state.p_cc);
    
    p_crit(2) = 0;
    m_dot_ox_crit(2) = 0;
    m_dot_lox(2) = 0;
    m_dot_gox(2) = m_dot_ox*...
        (state.m_gox)/(state.m_gox + state.m_oxtank_press);
    m_dot_oxtank_press(2) = m_dot_ox*...
        (state.m_oxtank_press)/(state.m_gox + state.m_oxtank_press);
    Q(2) = m_dot_ox*state.V_ox_ullage/(state.m_gox + state.m_oxtank_press);
    
    %% Total Flow Rate
    if state.m_lox > dm_lox_tol
        frac_lox = 1;
    else
        frac_lox = max(0,state.m_lox/dm_lox_tol);
    end
    
    m_dot_lox = frac_lox*m_dot_lox(1) + (1-frac_lox)*m_dot_lox(2);
    m_dot_gox = frac_lox*m_dot_gox(1) + (1-frac_lox)*m_dot_gox(2);
    m_dot_oxtank_press = frac_lox*m_dot_oxtank_press(1) + (1-frac_lox)*m_dot_oxtank_press(2);
    p_crit = frac_lox*p_crit(1) + (1-frac_lox)*p_crit(2);
    m_dot_ox_crit = frac_lox*m_dot_ox_crit(1) + (1-frac_lox)*m_dot_ox_crit(2);
    Q = frac_lox*Q(1) + (1-frac_lox)*Q(2);
else
    m_dot_lox = 0;
    m_dot_gox = 0;
    m_dot_oxtank_press = 0;
    p_crit = state.p_oxtank;
    m_dot_ox_crit = 0;
    Q = 0;
end


%% Oxidizer tank draining
dW = -(state.p_gox+state.p_oxtank_press)*Q; % work done on gas by liquid
T_dot_drain = dW/state.oxtank_m_cv;
end

function [m_dot_ox, m_dot_crit, p_crit] = LN2OMDot(inputs, A_inj, Pup, Tup, Pdown)
% Calculation of the mass flow rate of liquid nitrous oxide

% Set pressure to be a minimum of vapor pressure to address modeling
% limitation
n2o_prop = N2O_Properties(Tup);
if Pup < n2o_prop.Pvap
    Pup = n2o_prop.Pvap;
end

% Use 2-phase model for low supercharge pressure ratio
if Pup/n2o_prop.Pvap < 3
    [G, G_crit, p_down_crit] = FindG(Tup, Pup, Pdown);
else
    G = sqrt(2*(Pup-Pdown)*n2o_prop.rho_l);
    G_crit = 0;
    p_down_crit = 0;
end
m_dot_ox = inputs.ox.Cd_injector*A_inj*G;
m_dot_crit = inputs.ox.Cd_injector*A_inj*G_crit;
p_crit = p_down_crit;
end

function [m_dot_vap, T_dot_vap] = N2OTankEquilibrium(inputs, state)
%N2OTankEquilibrium Find equilibrium state of nitrous oxide tank
%   Inputs:
%       state: object describing motor characteristics
%       state: object representing state of motor
%   Outpus:
%       m_dot_vap: mass of nitrous oxide vaporizing from liquid to gas
%       T_dot_vap: time rate of change of temperature due to oxidizer
%           vaporization

% Constants
vap_const = 1e2; % vaporization constant
dpvap_dT = 5e4; % approximate slope of vapor pressure w.r.t. temperature

% Establish vapor equilibrium for oxidizer
m_dot_vap = vap_const*(state.N2O_properties.Pvap - state.p_gox)*...
    state.oxtank_m_cv/(state.N2O_properties.deltaE_vap*dpvap_dT); 
    % mass flow rate vaporizing (liquid to gas)
if m_dot_vap > 0 && state.m_lox < 0
    m_dot_vap = 0; 
        % no vaporization if liquid oxidizer is gone
elseif m_dot_vap < 0 && state.m_gox < 0
    m_dot_vap =  0;
end
H_dot_tank_vap = m_dot_vap*-state.N2O_properties.deltaH_vap;
T_dot_vap = H_dot_tank_vap/state.oxtank_m_cp;
end


function [F_thrust, OF, M_e, p_exit, p_shock, m_cc_dot, M_cc_dot, ...
    gamma_cc_dot, T_cc_dot] = CombustionChamber(inputs, state, m_dot_ox, m_dot_f)
%CombustionChamber Model combustion chamber dynamics
%   Inputs:
%       state: object describing motor characteristics
%       state: object representing state of system
%       m_dot_ox: oxidizer mass flow rate
%       m_dot_f: fuel mass flow rate
%   Outputs:
%       F_thrust: thrust force
%       OF: oxidizer to fuel mass ratio
%       M_e: nozzle exit mach number
%       p_exit: nozzle exit pressure
%       p_shock: nozzle exit pressure at which shock occurs
%       m_cc_dot: time rate of change of combustion chamber gas mass
%       M_cc_dot: time rate of change of combustion chamber gas molecular mass
%       gamma_cc_dot: time rate of change of combustion chamber gas ratio 
%           of specific heats
%       T_cc_dot: time time rate of change of combustion chamber gas
%           temperature


%% Combustion Dynamics
if m_dot_ox > 0
    OF = m_dot_ox/m_dot_f;
    [gamma_ex, iMW, R_ex, iTc, c_star] = ...
        CombustionCalc(OF,state.p_cc,inputs,inputs.comb_data);
else % If combustion chamber pressure is greater than tank pressure
    gamma_ex = 0;
    iMW = 0;
    iTc = 0;
    OF = 0;
    c_star = 1;
end               

%% Solve Flow Through Nozzle    
%Nozzle exit area calculation based on throat area, expansion ratio
area_throat = pi*(inputs.d_throat)^2/4;
exit_area = area_throat*inputs.exp_ratio;
d_exit = inputs.d_throat*sqrt(inputs.exp_ratio);

% Run nozzle calculations
if state.p_cc > inputs.p_amb
    [ ~, p_exit, p_shock, u_exit, ~, M_e ] = NozzleCalc( inputs.d_throat, ...
        d_exit, state.T_cc, state.p_cc, state.gamma_cc, state.M_cc, inputs.p_amb); 
    m_dot_ex = state.p_cc*area_throat/c_star;
else
    p_exit = inputs.p_amb;
    p_shock = 0;
    u_exit = 0;
    m_dot_ex = 0;
    M_e = 0;
end

% Account for nozzle efficiency
u_exit = u_exit*sqrt(inputs.nozzle_efficiency);

%Calculate thrust force
C_T = inputs.nozzle_correction_factor*(m_dot_ex*u_exit + ...
    (p_exit-inputs.p_amb)*exit_area)/(area_throat*state.p_cc);
F_thrust = c_star*C_T*m_dot_ex;

accel_net = F_thrust/(inputs.mass_dry_rocket + state.m_fuel + ...
    state.m_lox + state.m_gox);

%% Combustion chamber properties
m_cc_dot = m_dot_ox + m_dot_f - m_dot_ex;
M_cc_dot = (m_dot_ox + m_dot_f)*(iMW-state.M_cc)/state.m_cc;
gamma_cc_dot = (m_dot_ox + m_dot_f)*...
    (gamma_ex-state.gamma_cc)/state.m_cc;
T_cc_dot = (m_dot_ox + m_dot_f)*(iTc-state.T_cc)/state.m_cc;
end

function [m_dot_f, r_dot] = HybridMdotFuel(inputs, state, m_dot_ox)
%HybridMdotFuel Calculate mass flow rate of fuel and regression rate
%   Inputs:
%       inputs: object describing motor characteristics
%       state: object representing system state
%   Outputs:
%       m_dot_f: mass flow rate of fuel
%       r_dot: regression rate of fuel grain

areap = pi*state.rad_port^2;
r_rate = inputs.fuel.a*(m_dot_ox/areap)^(inputs.fuel.n);
r_vol = 2*pi*state.rad_port*r_rate*inputs.fuel.grain_length;
m_dot_f = r_vol * inputs.fuel.rho;
r_dot = r_rate;
end

function [m_dot_f, T_dot_drain] = FuelTankMDot(inputs, state, time)
%FuelTankMDot Simulates dynamics of draining fuel tank
%   Inputs:
%       - inputs: object representing motor inputs
%       - state: object representing state of system
%       - time: time
%   Outputs:
%       - m_dot_f: mass flow rate out of tank of liquid fuel
%       - T_dot_drain: fuel tank temperature change from draining

% If tank pressure is greater than combustion chamber pressure
injector_area = inputs.fuel.Cd_injector*inputs.fuel.injector_area*inputs.Throttle(time);
if (state.p_fueltank > state.p_cc)&&...
        (injector_area>0)
    G = sqrt(2*(state.p_fueltank-state.p_cc)*inputs.fuel.rho);
    m_dot_f = injector_area*G;
else
    m_dot_f = 0;
end
Q = m_dot_f/inputs.fuel.rho;

%% Oxidizer tank draining
dW = -(state.p_fueltank)*Q; % work done on gas by liquid
T_dot_drain = dW/(state.m_fueltank_press*inputs.fuel_pressurant.gas_properties.c_v);
end

function [gamma_ex, iMW, R_ex, iTc, c_star] = ...
    CombustionCalc(OF,Pc,inputs,CombData)
%CombustionCalc Calculate combustion properties. 
%   Inputs:
%       OF: oxidizer to fuel mass flow ratio
%       Pc: combustion chamber static pressure
%       inputs: motor input structure
%       CombDataFit: fit object with combustion properties over range of OF
%           and chamber pressure
%   Outputs:
%       gamma_ex: ratio of specific heats of combustion gases
%       iMW: molecular mass of combustion gases
%       R_ex: exhaust specific gas constant
%       iTc: stagnation temperature of combustion gases
%       c_star: characteristic velocity of combustion

%% Constants
R_u = 8.3144621; % Universal gas constant [J/mol*K]

%% Calculate Combustion Properties
% Correct OF ratio or chamber pressure if over boundaries
if OF < CombData.OF_range(1)
    OF = CombData.OF_range(1);
elseif OF > CombData.OF_range(2)
    OF = CombData.OF_range(2);
end
if Pc < CombData.Pc_range(1)
    Pc = CombData.Pc_range(1);
elseif Pc > CombData.Pc_range(2)
    Pc = CombData.Pc_range(2);
end

% Use empirical formula to calculate chamber temperature, exhaust
% properties
iTc = FastInterp2(CombData.OF, CombData.Pc, CombData.Tc, OF, Pc);
iMW = FastInterp2(CombData.OF, CombData.Pc, CombData.M, OF, Pc);
R_ex = R_u/iMW;
gamma_ex = FastInterp2(CombData.OF, CombData.Pc, CombData.gamma, OF, Pc);
c_star = FastInterp2(CombData.OF, CombData.Pc, CombData.c_star, OF, Pc);

% Apply c-star efficiency
c_star = c_star*inputs.c_star_efficiency;
iTc = iTc*inputs.c_star_efficiency^2;

% Check Tc, MW, gamma_ex for invalid values
if (iTc < 0) || ~isreal(iTc) || isnan(iTc)
    error('Invalid Tc calculated: Tc = %.3f', iTc);
end
if (iMW < 0) || ~isreal(iMW) || isnan(iMW)
    error('Invalid MW calculated: MW = %.3f', iMW);
end
if (gamma_ex < 0) || ~isreal(gamma_ex) || isnan(gamma_ex)
    error('Invalid gamma_ex calculated: gamma_ex = %.3f', gamma_ex);
end
end

function [T_dot_press, m_dot_press] = PressurantFlow(p_presstank, ...
    T_presstank, p_downtank, T_downtank, m_cv_downtank, pressurant)
%PressurantFlow Models flow of pressurant gas from source tank to
%propellant tank
%   Inputs:
%       p_presstank: pressure in source tank
%       T_presstank: temperature in source tank
%       p_downtank: pressure in propellant tank
%       T_downtank: temperature in propellant tank
%       m_cv_downtank: thermal capacity (constant volume) of propellant 
%           tank
%       pressurant: object with properties of pressurant gas
%   Outputs:
%       T_dot_press: temperature time rate of change of propellant tank 
%           due to pressurant inflow
%       m_dot_press: mass flow rate of pressurant gas from source tank to
%           propellant tank

if (p_presstank > pressurant.set_pressure) && ...
        (p_downtank < pressurant.set_pressure) && (pressurant.tank_volume > 0)
    % Determine mass flow of pressurant
    d_reg = sqrt(4/pi*(pressurant.flow_CdA));
    [ ~, ~, ~, ~, m_dot_press, ~ ] = NozzleCalc( d_reg, d_reg, ...
        T_presstank, p_presstank, ...
        pressurant.gas_properties.gamma, ...
        pressurant.gas_properties.molecular_mass, ...
        p_downtank);
    % Taper flow rate to zero within pressure tolerance
    dp_tol = 1e5;  % pressure difference tolerance
    m_dot_press = min(m_dot_press,m_dot_press*...
        (pressurant.set_pressure - p_downtank)/dp_tol);
    
    % Update tank conditions based on more pressurant flow
    T_press_down = T_presstank*(p_downtank/p_presstank)^...
        ((pressurant.gas_properties.gamma-1)/pressurant.gas_properties.gamma);
    rho_press_tank = p_downtank/(pressurant.gas_properties.R_specific*T_press_down);
    dW = p_downtank*m_dot_press/rho_press_tank;
    
    % Thermal equilibrium within tank
    T_dot_press = (dW + (T_press_down-T_downtank)*...
        m_dot_press*pressurant.gas_properties.c_v)/m_cv_downtank;
else
    m_dot_press = 0;
    T_dot_press = 0;
end
end

function p = PVDW(T, rho, R, a, b)
%PVDW VanDerWaal's equation of state for pressure
    p = R*T./(1./rho-b)-a*rho.^2;
end

function T = TVDW(p, rho, R, a, b)
%PVDW VanDerWaal's equation of state for temperature
    T = (p + a*rho).*(1./rho - b)/R;
end

function rho = RhoVDW(p, T, R, a, b)
%PVDW VanDerWaal's equation of state for density, gas values
    vdw_roots = roots([-a*b, a, - p*b - R*T, p]);
    real_vdw_roots = [];
    for ii = 1:length(vdw_roots)
        if isreal(vdw_roots(ii))
            real_vdw_roots(end+1) = vdw_roots(ii);
        end
    end
    rho = min(real_vdw_roots);
end
