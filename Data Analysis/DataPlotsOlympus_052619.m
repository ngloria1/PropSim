%% Post-Processing Script
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

%% Import data and convert to SI units
data_fuel = readtable("../Test Data/2019_05_26_Coldflow_A.csv");
data_ox = readtable("../Test Data/2019_05_26_Coldflow_B.csv");

% Time offsets
fuel_start_time = 225.7;
ox_start_time = 996.4;
test_duration = 58;

time_fuel = data_fuel.Time*1e-3 - fuel_start_time;
time_ox = data_ox.Time*1e-3 - ox_start_time;
pft = data_fuel.FuelTank*psi_to_Pa; % Pressure Fuel Tank
pot = data_ox.x4_OxTank*psi_to_Pa; % Pressure Ox Tank
pos = data_ox.x1_NitrousSupply*psi_to_Pa; % Pressure Ox Supply
pom = data_ox.x5_Manifold*psi_to_Pa; % Pressure Ox Supply

time = union(time_fuel,time_ox);
unique_time_fuel = diff(time_fuel) ~= 0;
unique_time_ox = diff(time_ox) ~= 0;
pft = interp1(time_fuel(unique_time_fuel),pft(unique_time_fuel),time);
pot = interp1(time_ox(unique_time_ox),pot(unique_time_ox),time);
pos = interp1(time_ox(unique_time_ox),pos(unique_time_ox),time);
pom = interp1(time_ox(unique_time_ox),pom(unique_time_ox),time);
pft_sliced = pft(time > 0 & time < test_duration);
pot_sliced = pot(time > 0 & time < test_duration);
pom_sliced = pom(time > 0 & time < test_duration);
test_time = time(time > 0 & time < test_duration);

%% User Input
rho_fuel = 1000; % kg/m^3
t_fuel_cutoff = 38.4;
p_amb = 14.7*psi_to_Pa; % ambient pressure
gamma_fuelpress = 1.4; % N2
V_ullage_fuel = 0.88e-3; % m^3
pft_start = 603*psi_to_Pa; % psi, after valve open
pft_shutoff = 272*psi_to_Pa; % psi
dm_fuel = 1e-3*rho_fuel; % fuel mass expended, kg

%% Do Calculations
dp_fuel = pft_sliced - p_amb;
tau_CdA = 1.0; % time span over which to measure CdA, s
pft_filt = SimpleFilter(time,pft,tau_CdA);
pft_filt = pft_filt(time > 0 & time < test_duration);
dp_fuel_filt = pft_filt - p_amb;
fuel_flow_time = test_time < t_fuel_cutoff;
% Calculate fuel mass flow based on ullage expansion
m_fuel = rho_fuel * V_ullage_fuel*((pft_start./pft_filt).^(1/gamma_fuelpress)-1);
m_dot_fuel = DiffFilter(test_time,m_fuel,tau_CdA);
fuel_scaling_factor = dm_fuel/trapz(test_time(fuel_flow_time),m_dot_fuel(fuel_flow_time));
m_dot_fuel = fuel_scaling_factor*m_dot_fuel;
m_fuel = fuel_scaling_factor*m_fuel;
% Calculate fuel CdA
fuel_CdA = m_dot_fuel./(sqrt(2*rho_fuel)*sqrt(dp_fuel_filt));
fuel_CdA_int = dm_fuel./(sqrt(2*rho_fuel)*...
    trapz(test_time(fuel_flow_time),sqrt(dp_fuel(fuel_flow_time))));

%% Plots
figure()
plot(time + fuel_start_time,pft/psi_to_Pa);
xlabel('Time [s]')
ylabel('Fuel Tank Pressure [psi]')
title('Fuel Tank Pressure Overview')

figure()
plot(time + ox_start_time,pot/psi_to_Pa,time + ox_start_time,pos/psi_to_Pa,...
    time + ox_start_time,pom/psi_to_Pa);
xlabel('Time [s]')
ylabel('Pressure [psi]')
title('Oxidizer Test Overview')
legend({'Rocket Oxidizer Tank','Supply Oxidizer','Oxidizer Manifold'})

figure()
plot(test_time,pot_sliced/psi_to_Pa,test_time,pom_sliced/psi_to_Pa,...
    test_time,pft_sliced/psi_to_Pa);
xlabel('Time [s]')
ylabel('Pressure [psi]')
title('Engine Pressures Overview')
legend({'Oxidizer Tank','Oxidizer Manifold','Fuel Tank'})

figure()
plot(test_time(fuel_flow_time),m_fuel(fuel_flow_time))
xlabel('Time (s)');
ylabel('Fuel Mass Consumed (kg)');

figure()
plot(test_time(fuel_flow_time),m_dot_fuel(fuel_flow_time))
xlabel('Time (s)');
ylabel('Fuel Mass Flow Rate (kg/s)');

figure()
plot(test_time(fuel_flow_time),fuel_CdA(fuel_flow_time)*1e6)
ylim([0 2*median(fuel_CdA(fuel_flow_time))*1e6])
xlabel('Time (s)');
ylabel('Fuel CdA (mm^2)');

fprintf('Calculations:\n')
fprintf('Fuel Mass Spent: %.3g kg (%.3g lbm)\n', dm_fuel, dm_fuel/lbm_to_kg)
fprintf('Fuel CdA: %.3g mm^2\n', fuel_CdA_int*1e6)

%% Export Data
ExportData(test_time,pft_sliced,pom_sliced,pot_sliced,...
    zeros(size(time)),zeros(size(time)),zeros(size(time)))

%% Subfunctions
function ExportData(test_time,pft,pom,pot,we,ft,pcc)
    [file,path] = uiputfile('*.mat');
    save([path file],'test_time','pft','pom','pot','we','ft','pcc');
end

function time_fixed = FixTimeOverflow(time)
%FixTimeOverflow Rectify overflow time indices

time_fixed = time;
dt = diff(time);
median_dt = median(dt);
overflow_indices = find(abs(diff(time)) > 10 * median_dt);
diff_overflow = dt(overflow_indices);
for ii = 1:length(overflow_indices)
    time_fixed((overflow_indices(ii)+1):end) = ...
        time_fixed((overflow_indices(ii)+1):end) - diff_overflow(ii) + median_dt;
end
end

function smoothed = SimpleFilter(time,values,tau)
dt = median(diff(time));
Wn = dt/(tau);
values(isnan(values)) = 0;
[b, a] = butter(1,Wn,'low');
smoothed = filter(b, a, values);
end

function smoothed = DiffFilter(time,values,tau)
dt = median(diff(time));
dn = tau/dt;
smoothed = (circshift(values,round(0.5*dn)) - circshift(values,-round(0.5*dn)))./...
     (circshift(time,round(0.5*dn)) - circshift(time,-round(0.5*dn)));
end