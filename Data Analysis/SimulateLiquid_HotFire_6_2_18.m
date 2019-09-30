%% Run Performance Code

cd('..')
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
test_data.test_plots_on = 1; % Import tests data and plot against simulation data
test_data.test_data_file = '6_02_18.mat'; % File from which to import test data
test_data.t_offset = 0; % Time offset of test data wrt simulation data [s]

% 1: simulate combustion (hot fire), 0: no combustion (cold flow)
mode.combustion_on = 1;
% 1: simulate flight conditions (i.e. acceleration head), 0: ground test 
% conditions
mode.flight_on = 0;
% 'hybrid' for solid fuel, liquid oxidizer, 'liquid' for liquid fuel and
% oxidizer
mode.type = 'liquid';

%% Input Parameters

inputs.CombustionData = fullfile('Combustion Data', 'CombustionData_T1_N2O.mat');

%-------Gases-----------------------

helium = Gas();
helium.c_v = 3.12e3; % J/kg*K
helium.molecular_mass = 4.0026e-3; % kg/mol

nitrogen = Gas();
nitrogen.c_v = 0.743e3; % J/kg*K
nitrogen.molecular_mass = 2*14.0067e-3; % kg/mol

%-------Injector Properties----------

%Injector Exit Area
inputs.ox.injector_area = 28.5*0.95*mm_to_m^2;
inputs.fuel.injector_area = 2.1*mm_to_m^2;

% Ball Valve Time to Injector Area (s)
inputs.dt_valve_open = 0.01;

%Discharge Coefficient
inputs.ox.Cd_injector = 1;
inputs.fuel.Cd_injector = 1;

%-------Rocket Properties--------
%Rocket Dry Mass
inputs.mass_dry_rocket = 50*lbm_to_kg;

%-------Oxidizer Properties--------
%Tank Volume
inputs.ox.V_tank = 6.63*L_to_m3; 

%Nitrous Volume
inputs.ox.V_l = 6.0*L_to_m3; 

%Tank Inner Diameter
inputs.ox.tank_id = 3.75*in_to_m;

%Distance from Bottom of Tank to Injector
inputs.ox.h_offset_tank = 0*in_to_m;

%Main Flow Line Diameter
inputs.ox.d_flowline = .25*in_to_m;

%Tank Temperature (K)
inputs.ox.T_tank = 290; %287.3

%-------Oxidizer Pressurant Properties--------

inputs.ox_pressurant = Pressurant('oxidizer');
inputs.ox_pressurant.gas_properties = helium;
inputs.ox_pressurant.set_pressure = 800*psi_to_Pa;
inputs.ox_pressurant.storage_initial_pressure = 4500*psi_to_Pa;
inputs.ox_pressurant.tank_volume = 3.5*L_to_m3;
inputs.ox_pressurant.flow_CdA = 8*mm_to_m^2;

%Are You Supercharging? (0 for 'NO' 1 for 'YES')
inputs.ox_pressurant.active = 0;

%-------Fuel Properties--------

%Tank Volume
inputs.fuel.V_tank = 1.88*L_to_m3; 

%Nitrous Volume
inputs.fuel.V_l = 1.00*L_to_m3; 

%Tank Inner Diameter
inputs.fuel.tank_id = 3.75*in_to_m;

%Distance from Bottom of Tank to Injector
inputs.fuel.h_offset_tank = 24*in_to_m;

%Main Flow Line Diameter(in)
inputs.fuel.d_flowline = .25*in_to_m;

inputs.fuel.rho = 795; %Kg/m^3

%-------Fuel Pressurant Properties--------

inputs.fuel_pressurant = Pressurant('fuel');
inputs.fuel_pressurant.gas_properties = nitrogen;
inputs.fuel_pressurant.set_pressure = 715*psi_to_Pa;
inputs.fuel_pressurant.storage_initial_pressure = 4500*psi_to_Pa;
inputs.fuel_pressurant.tank_volume = 0.0*L_to_m3;
inputs.fuel_pressurant.flow_CdA = 8*mm_to_m^2;

%Are You Supercharging? (0 for 'NO' 1 for 'YES')
inputs.fuel_pressurant.active = 1;

%-------Other Properties--------

%Combustion chamber dimensions
inputs.length_cc = 4*in_to_m;
inputs.d_cc = 3.75*in_to_m;

%Estimated nozzle efficiency
inputs.nozzle_efficiency = 0.875;
inputs.nozzle_correction_factor = 0.9830;

% Estimated combustion efficiency
inputs.c_star_efficiency = 0.95;

%Nozzle Throat diameter
inputs.d_throat = 2.388e-2;

%Expansion Ratio
inputs.exp_ratio = 3.5;

%Ambient Temperature
inputs.T_amb = 280;

%Ambient Pressure
inputs.p_amb = 12.5*psi_to_Pa;

% Load Combustion Data
inputs.comb_data = load(inputs.CombustionData); 
inputs.comb_data = inputs.comb_data.CombData;

%% Run Performance Code
PerformanceCode(inputs, mode, test_data);
cd('Data Analysis')
