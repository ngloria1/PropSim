%% Run Design Liquid
% This script will find design parameters to match desired characteristics
% (i.e. thrust, OF, Pc, etc.)

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

%% Goal/Design Parameters
goal.max_thrust = 400*lbf_to_N;
goal.OF = 8.0;
goal.total_impulse = 18e3;
goal.min_fuel_dp = 0.25; % min dp as % of tank pressure
goal.min_ox_dp = 0.35; % min dp as % of tank pressure
goal.ox_to_fuel_time = 1.00; % ratio of liquid oxidizer flow time to liquid fuel flow time

design.p_tanks = 650*psi_to_Pa;
design.ox_ullage = 0.1;
design.exp_ratio = 3.5;

%% Initial Input Parameters

initial_inputs.CombustionData = fullfile('Combustion Data', 'CombustionData_T1_N2O.mat');

%-------Gases-----------------------

helium = Gas();
helium.c_v = 3.12e3; % J/kg*K
helium.molecular_mass = 4.0026e-3; % kg/mol

nitrogen = Gas();
nitrogen.c_v = 0.743e3; % J/kg*K
nitrogen.molecular_mass = 2*14.0067e-3; % kg/mol

%-------Injector Properties----------

%Injector Exit Area
initial_inputs.ox.injector_area = 31.2*mm_to_m^2;
initial_inputs.fuel.injector_area = 2.31*mm_to_m^2;

% Ball Valve Time to Injector Area (s)
initial_inputs.dt_valve_open = 0.1;

%Discharge Coefficient
initial_inputs.ox.Cd_injector = 1;
initial_inputs.fuel.Cd_injector = 1;

%-------Rocket Properties--------
%Rocket Dry Mass
initial_inputs.mass_dry_rocket = 50*lbm_to_kg;

%-------Oxidizer Properties--------
%Tank Volume
initial_inputs.ox.V_tank = 12.7*L_to_m3; 

%Nitrous Volume
initial_inputs.ox.V_l = 11.4*L_to_m3; 

%Tank Inner Diameter
initial_inputs.ox.tank_id = 4.75*in_to_m;

%Distance from Bottom of Tank to Injector
initial_inputs.ox.h_offset_tank = 0*in_to_m;

%Main Flow Line Diameter
initial_inputs.ox.d_flowline = .25*in_to_m;

%Tank Temperature (K)
initial_inputs.ox.T_tank = 288;

%-------Oxidizer Pressurant Properties--------

initial_inputs.ox_pressurant = Pressurant('oxidizer');
initial_inputs.ox_pressurant.gas_properties = helium;
initial_inputs.ox_pressurant.set_pressure = 800*psi_to_Pa;
initial_inputs.ox_pressurant.storage_initial_pressure = 4500*psi_to_Pa;
initial_inputs.ox_pressurant.tank_volume = 3.5*L_to_m3;
initial_inputs.ox_pressurant.flow_CdA = 8*mm_to_m^2;

%Are You Supercharging? (0 for 'NO' 1 for 'YES')
initial_inputs.ox_pressurant.active = 0;

%-------Fuel Properties--------

%Tank Volume
initial_inputs.fuel.V_tank = 3.79*L_to_m3; 

%Nitrous Volume
initial_inputs.fuel.V_l = 1.39*L_to_m3; 

%Tank Inner Diameter
initial_inputs.fuel.tank_id = 4.75*in_to_m;

%Distance from Bottom of Tank to Injector
initial_inputs.fuel.h_offset_tank = 24*in_to_m;

%Main Flow Line Diameter(in)
initial_inputs.fuel.d_flowline = .25*in_to_m;

initial_inputs.fuel.rho = 786; %Kg/m^3

%-------Fuel Pressurant Properties--------

initial_inputs.fuel_pressurant = Pressurant('fuel');
initial_inputs.fuel_pressurant.gas_properties = nitrogen;
initial_inputs.fuel_pressurant.set_pressure = 625*psi_to_Pa;
initial_inputs.fuel_pressurant.storage_initial_pressure = 4500*psi_to_Pa;
initial_inputs.fuel_pressurant.tank_volume = 0.0*L_to_m3;
initial_inputs.fuel_pressurant.flow_CdA = 8*mm_to_m^2;

%Are You Supercharging? (0 for 'NO' 1 for 'YES')
initial_inputs.fuel_pressurant.active = 1;

%-------Other Properties--------

%Combustion chamber dimensions
initial_inputs.length_cc = 8*in_to_m;
initial_inputs.d_cc = 3.75*in_to_m;

% Estimated nozzle efficiency
initial_inputs.nozzle_efficiency = 0.95;
initial_inputs.nozzle_correction_factor = 0.9830;

% Estimated combustion efficiency
initial_inputs.c_star_efficiency = 0.85;

% Nozzle Throat diameter
initial_inputs.d_throat = 2.545e-2;

% Ambient Temperature
initial_inputs.T_amb = 280;

% Ambient Pressure
initial_inputs.p_amb = 12.74*psi_to_Pa;

% Load Combustion Data
initial_inputs.comb_data = load(initial_inputs.CombustionData); 
initial_inputs.comb_data = initial_inputs.comb_data.CombData;


%% Run Code
DesignLiquid(initial_inputs, goal, design, true);
