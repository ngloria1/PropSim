%% Run Performance Code

addpath(fullfile(pwd, 'Supporting Functions'))

clear
close all

% Unit Conversion
psi_to_Pa = 6894.75729; % 1 psi in Pa
in_to_m = 0.0254; % 1 in in m
mm_to_m = 1e-3; % 1 mm in m
lbf_to_N = 4.44822162; % 1 lbf in N
lbm_to_kg = 0.453592; % 1 lbm in kg
atm_to_Pa = 101325; % 1 atm in Pa
L_to_m3 = 1e-3; % 1 L in m^3


%% Options
test_data.test_plots_on = 0; % Import tests data and plot against simulation data
test_data.test_data_file = '5-28-16.mat'; % File from which to import test data
test_data.t_offset = 0; % Time offset of test data wrt simulation data [s]

% 1: simulate combustion (hot fire), 0: no combustion (cold flow)
mode.combustion_on = 1;
% 1: simulate flight conditions (i.e. acceleration head), 0: ground test 
% conditions
mode.flight_on = 0;
% 'hybrid' for solid fuel, liquid oxidizer, 'liquid' for liquid fuel and
% oxidizer
mode.type = 'hybrid';

%% Input Parameters

inputs.CombustionData = fullfile('Combustion Data', 'CombustionData_80PR_20PP.mat');

%-------Injector Properties----------

%Injector Exit Area
inputs.ox.injector_area = 64.6*mm_to_m^2;

% Ball Valve Time to Injector Area (s)
inputs.dt_valve_open = 0.1;

%Discharge Coefficient
inputs.ox.Cd_injector = 1;

%-------Rocket Properties--------
%Rocket Dry Mass
inputs.mass_dry_rocket = 70*lbm_to_kg;

%Tank Volume
inputs.ox.V_tank = 20*L_to_m3; 

%Nitrous Volume
inputs.ox.V_l = 18*L_to_m3; 

%Tank Inner Diameter
inputs.ox.tank_id = 5.375*in_to_m;

%Distance from Bottom of Tank to Injector
inputs.ox.h_offset_tank = 6*in_to_m;

%Main Flow Line Diameter(in)
inputs.ox.d_flowline = .5*in_to_m;

%Tank Temperature (K)
inputs.ox.T_tank = 289.81;

%-------Hybrid Properties--------

%Fuel Density
inputs.fuel.rho = 880; %Kg/m^3

%Grain Length
inputs.fuel.grain_length = 16.5*in_to_m;

%Grain Diameter
inputs.fuel.grain_od = (5-0.160)*in_to_m;

%Grain geometry
inputs.fuel.port_rad = 1.2*in_to_m;

%Regression Rate Coefficients
inputs.fuel.n = 0.67; % []
inputs.fuel.a = 0.5*0.104*mm_to_m;

%-------Other Properties--------

%Combustion chamber dimensions
inputs.length_cc = 26*in_to_m;
inputs.d_cc = 5.375*in_to_m;

%Estimated nozzle efficiency
inputs.nozzle_efficiency = 0.95;
inputs.nozzle_correction_factor = 0.9830;

% Estimated combustion efficiency
inputs.c_star_efficiency = 0.87;

%Nozzle Throat diameter
inputs.d_throat = 3.77e-2;

%Expansion Ratio
inputs.exp_ratio = 3.9;

%Ambient Temperature
inputs.T_amb = 280;

%Ambient Pressure
inputs.p_amb = 12.5*psi_to_Pa;

%-------Pressurant Properties--------

helium = Gas();
helium.c_v = 3.12; % J/kg*K
helium.molecular_mass = 4.0026e-3; % kg/mol

inputs.ox_pressurant = Pressurant('oxidizer');
inputs.ox_pressurant.gas_properties = helium;
inputs.ox_pressurant.set_pressure = 800*psi_to_Pa;
inputs.ox_pressurant.storage_initial_pressure = 4500*psi_to_Pa;
inputs.ox_pressurant.tank_volume = 3.5*L_to_m3;
inputs.ox_pressurant.flow_CdA = 8e2*(1e-3)^3;

%Are You Supercharging? (0 for 'NO' 1 for 'YES')
inputs.ox_pressurant.active = 0;

inputs.comb_data = load(inputs.CombustionData); 
inputs.comb_data = inputs.comb_data.CombData;

%% Run Performance Code
PerformanceCode(inputs, mode, test_data);
