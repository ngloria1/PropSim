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

data = readtable("../Test Data/2019_06_08_Hotfire.csv");

%% User Input
t_shutoff_fuel = 2.25; % s, time at shutoff of fuel
rho_fuel = 795; % kg/m^3
gamma_fuelpress = 1.4; % N2
V_ullage_fuel = 1.465e-3; % m^3
pft_start = 458*psi_to_Pa; % psi, after valve open
pft_shutoff = 350*psi_to_Pa; % psi
A_star = pi/4*(2.545e-2)^2; % m^2
m_ox = 3.30; % oxidizer mass expended, kg (3.30 overall, 2.94 liquid)
dm_fuel = 0.415e-3*rho_fuel; % fuel mass expended, kg

c_star_theo = 1506; % m/s
C_f_theo = 1.38; 
Isp_theo = c_star_theo*C_f_theo; % m/s

%% Import data and convert to SI units
time = data.Time*1e-3;
pom = data.x5_Manifold*psi_to_Pa;
pft = data.FuelTank*psi_to_Pa;
pot = data.x4_OxTank*psi_to_Pa;
ft = -data.loadCell*lbf_to_N* 3.3 * 1000 / (1024 * 50.1 * 0.015);
pcc = data.x2_NitrousHeatXger*psi_to_Pa*3/4;
dp_fuel = pft - pcc;
dp_ox = pom - pcc;

%% Find burn end and start times
burn_start_time = 1220.45;
burn_end_time = 1220.45 + 7.8;
burn_time = 3.6;

%% Generate Selection Indices
test_start_time = burn_start_time - 5;
test_end_time = burn_end_time + 30;
burn_start_index = find(time>burn_start_time,1, 'first');
burn_end_index = find(time>burn_end_time, 1, 'first');
test_start_index = find(time>test_start_time,1, 'first');
test_end_index = find(time>test_end_time, 1, 'first');

burn_ind = burn_start_index:burn_end_index;
test_ind = test_start_index:test_end_index;

% Mark T0 as burn start
time = time - burn_start_time;

% Time-based indices
nominal_flow_ind = time > 0 & time < t_shutoff_fuel;

%% Tare load cell
ft = ft - mean(ft(time < 0 & time > -5));

%% Filtering
tau_filter = 1.0; % s
pft_filt = SimpleFilter(time,pft,tau_filter);
dp_fuel_filt = SimpleFilter(time,dp_fuel,tau_filter);

%% Do Calculations
% Calculate fuel mass flow based on ullage expansion
m_fuel = rho_fuel * V_ullage_fuel*((pft_start./pft).^(1/gamma_fuelpress)-1);
dt = median(diff(time(nominal_flow_ind)));
tau_CdA = 1.0; % time span over which to measure CdA, s
m_dot_fuel = DiffFilter(time,m_fuel,tau_CdA);
dm_fuel_test = rho_fuel * V_ullage_fuel*((pft_start/pft_shutoff)^(1/gamma_fuelpress)-1);
fuel_scaling_factor = dm_fuel/(rho_fuel * V_ullage_fuel*((pft_start/pft_shutoff)^(1/gamma_fuelpress)-1));
m_dot_fuel = fuel_scaling_factor*m_dot_fuel;
% Calculate fuel CdA
fuel_CdA = m_dot_fuel./(sqrt(2*rho_fuel)*sqrt(dp_fuel_filt));
fuel_CdA_int = dm_fuel./(sqrt(2*rho_fuel)*...
    trapz(time(nominal_flow_ind),sqrt(dp_fuel(nominal_flow_ind))));

% Calculate Performance Characteristics
Impulse = trapz(time(burn_ind),ft(burn_ind));
pcc_int = trapz(time(burn_ind),pcc(burn_ind));
C_f = ft./(pcc*A_star);
C_f_int = Impulse/(pcc_int*A_star);
c_star_int = (pcc_int*A_star)/(m_ox + dm_fuel_test);
Isp = Impulse/(m_ox + dm_fuel_test);
OF = m_ox/dm_fuel_test;
ft_avg = Impulse/burn_time;
pcc_avg = pcc_int/burn_time;

fprintf('Calculations:\n')
fprintf('Impulse: %.3g kN*s (%.3g lbf*s)\n', Impulse/1e3, Impulse/lbf_to_N)
fprintf('Burn Time: %.3g s\n', burn_time)
fprintf('Avg. Thrust: %.3g kN (%.3g lbf)\n', ft_avg/1e3, ft_avg/lbf_to_N)
fprintf('Avg. Pcc: %.3g MPa (%.3g psi)\n', pcc_avg/1e6, pcc_avg/psi_to_Pa)
fprintf('Isp: %.3g m/s (%.3g s)\n', Isp, Isp/9.81)
fprintf('Isp efficiency: %.3g %%\n', Isp/Isp_theo)
fprintf('C_f: %.3g \n', C_f_int)
fprintf('C_f efficiency: %.3g %%\n', C_f_int/C_f_theo)
fprintf('C*: %.3g m/s (%.3g s)\n', c_star_int, c_star_int/9.81)
fprintf('C* efficiency: %.3g %%\n', c_star_int/c_star_theo)
fprintf('OF (avg): %.3g\n', OF)
fprintf('Fuel Mass Spent: %.3g kg (%.3g lbm)\n', dm_fuel, dm_fuel/lbm_to_kg)
fprintf('Ox. Mass Spent: %.3g kg (%.3g lbm)\n', m_ox, m_ox/lbm_to_kg)
fprintf('Fuel CdA: %.3g mm^2\n', fuel_CdA_int*1e6)

%% Generate Plots
% Test Time Plots
figure()
plot((time(test_ind)), pot(test_ind)/psi_to_Pa,...
    (time(test_ind)), pom(test_ind)/psi_to_Pa,...
    (time(test_ind)), pft(test_ind)/psi_to_Pa,...
    (time(test_ind)), pcc(test_ind)/psi_to_Pa)
xlabel('Time (s)');
ylabel('Pressures (psi)');
legend({'Ox. Tank','Ox. Manifold','Fuel Tank','Combustion Chamber'})
figure()
plot((time(time > 0 & time < t_shutoff_fuel)), ...
    1e2*dp_fuel(time > 0 & time < t_shutoff_fuel)./pft(time > 0 & time < t_shutoff_fuel),...
    (time(time > 0 & time < t_shutoff_fuel)), ...
    1e2*dp_ox(time > 0 & time < t_shutoff_fuel)./pot(time > 0 & time < t_shutoff_fuel))
xlabel('Time (s)');
ylabel('Pressure Drop (% of tank pressure)');
legend({'Fuel','Oxidizer'})
figure()
plot((time(test_ind)), ft(test_ind)/lbf_to_N)
xlabel('Time (s)');
ylabel('Thrust (lbf)');

% Performance Parameters
figure()
plot(time(burn_ind), C_f(burn_ind))
xlabel('Time (s)');
ylabel('Thrust Coefficient ()');

% CdA Plots
figure()
plot(time(time > -1 & time < t_shutoff_fuel + 10),pft_filt(time > -1 & time < t_shutoff_fuel + 10)/psi_to_Pa,...
    time(time > -1 & time < t_shutoff_fuel + 10),pft(time > -1 & time < t_shutoff_fuel + 10)/psi_to_Pa)
xlabel('Time (s)');
ylabel('Pressure Fuel Tank - filtered (psi)');
figure()
plot(time(time > -1 & time < t_shutoff_fuel),m_fuel(time > -1 & time < t_shutoff_fuel))
xlabel('Time (s)');
ylabel('Fuel Mass Consumed (kg)');
figure()
plot(time(time > -1 & time < t_shutoff_fuel),m_dot_fuel(time > -1 & time < t_shutoff_fuel))
xlabel('Time (s)');
ylabel('Fuel Mass Flow Rate (kg/s)');
figure()
plot(time(nominal_flow_ind),fuel_CdA(nominal_flow_ind)*1e6)
ylim([0 2*median(fuel_CdA(nominal_flow_ind))*1e6])
xlabel('Time (s)');
ylabel('Fuel CdA (mm^2)');

%% Export Data
ExportData(time(test_ind),pft(test_ind),pom(test_ind),pot(test_ind),...
    zeros(size(time(test_ind))),ft(test_ind),pcc(test_ind))

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