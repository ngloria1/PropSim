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

%% Calibration

data = readtable("../Test Data/2019_06_17_Calibration.csv"); 

time = data.Channel8*1e-3;
pom = data.Channel3*psi_to_Pa;
pft = data.Channel4*psi_to_Pa;
pot = data.Channel5*psi_to_Pa;
ft = -data.Channel7*lbf_to_N* 3.3 * 1000 / (1024 * 50.1 * 0.015);
pcc = data.Channel6*psi_to_Pa;
dp_fuel = pft - pcc;
dp_ox = pom - pcc;

ft_filt = SimpleFilter(time,ft,5.0,'low');

figure
plot(time,pcc/psi_to_Pa,time,pom/psi_to_Pa)
legend({'Combustion Chamber', 'Ox. Inj. Manifold'})
xlabel('Time [s]')
ylabel('Pressure [psi]')
title('Pressurized to 800 psi')

pcc_factor = 800/713;
pom_factor = 800/713;

figure
plot(time,ft_filt/lbf_to_N+3347)
xlabel('Time [s]')
ylabel('Thrust Load Cell [lbf]')

time_unloaded = [480, 500];
time_loaded = [520, 550];
ft_unloaded = mean(ft(time > time_unloaded(1) & time < time_unloaded(2)));
ft_loaded = mean(ft(time > time_loaded(1) & time < time_loaded(2)));
ft_factor = 14.45*9.81/(ft_loaded-ft_unloaded);

%% Test data

data = readtable("../Test Data/2019_06_14_Hotfire.csv");

%% User Input
t_shutoff_fuel = 10.05; % s, time at shutoff of fuel
rho_fuel = 795; % kg/m^3
gamma_fuelpress = 1.4; % N2
V_ullage_fuel = 2.4e-3; % m^3
pft_start = 658*psi_to_Pa; % psi, after valve open
pft_shutoff = 370*psi_to_Pa; % psi
A_star = pi/4*(2.545e-2)^2; % m^2
m_ox = 6.60; % oxidizer mass expended, kg (3.30 overall, 2.94 liquid)
dm_fuel = 1.42e-3*rho_fuel; % fuel mass expended, kg

c_star_theo = 1574; % m/s
C_f_theo = 1.4085; 
Isp_theo = c_star_theo*C_f_theo; % m/s

%% Import data and convert to SI units
time = data.Var8*1e-3;
pom = data.Var3*psi_to_Pa;
pft = data.Var4*psi_to_Pa;
pot = data.Var5*psi_to_Pa;
ft = data.Var7*lbf_to_N* 3.3 * 1000 / (1024 * 50.1 * 0.015)*ft_factor;
pcc = data.Var6*psi_to_Pa;
dp_fuel = pft - pcc;
dp_ox = pom - pcc;
pcc_diff = SimpleFilter(time,pcc,1.0,'high');

%% Find burn end and start times
burn_start_time = 2032.8;
burn_end_time = 2032.8 + 10.1;
burn_time = burn_end_time - burn_start_time;

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
pft_filt = SimpleFilter(time,pft,tau_filter,'low');
dp_fuel_filt = SimpleFilter(time,dp_fuel,tau_filter,'low');

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
pcc_rms_int = sqrt(trapz(time(burn_ind),pcc_diff(burn_ind).^2));
C_f = ft./(pcc*A_star);
C_f_int = Impulse/(pcc_int*A_star);
ft_factor = 0.95*C_f_theo/C_f_int;
ft = ft*ft_factor;
Impulse = Impulse*ft_factor;
C_f = C_f*ft_factor;
C_f_int = C_f_int*ft_factor;

c_star_int = (pcc_int*A_star)/(m_ox + dm_fuel_test);
Isp = Impulse/(m_ox + dm_fuel_test);
OF = m_ox/dm_fuel_test;
ft_avg = Impulse/burn_time;
pcc_avg = pcc_int/burn_time;

fprintf('Calculations:\n')
fprintf('Impulse: %.1f kN*s (%.0f lbf*s)\n', Impulse/1e3, Impulse/lbf_to_N)
fprintf('Burn Time: %.1f s\n', burn_time)
fprintf('Avg. Thrust: %.2f kN (%.0f lbf)\n', ft_avg/1e3, ft_avg/lbf_to_N)
fprintf('Avg. Pcc: %.2f MPa (%.0f psi)\n', pcc_avg/1e6, pcc_avg/psi_to_Pa)
fprintf('Avg. Pcc RMS: %.2f %%\n', pcc_rms_int/pcc_int*1e2)
fprintf('Isp: %.0f m/s (%.0f s)\n', Isp, Isp/9.81)
fprintf('Isp efficiency: %.1f %%\n', Isp/Isp_theo*1e2)
fprintf('C_f: %.2f \n', C_f_int)
fprintf('C_f efficiency: %.1f %%\n', C_f_int/C_f_theo*1e2)
fprintf('C*: %.0f m/s (%.1f s)\n', c_star_int, c_star_int/9.81)
fprintf('C* efficiency: %.1f %%\n', c_star_int/c_star_theo*1e2)
fprintf('OF (avg): %.1f\n', OF)
fprintf('Fuel Mass Spent: %.2f kg (%.2f lbm)\n', dm_fuel, dm_fuel/lbm_to_kg)
fprintf('Ox. Mass Spent: %.2f kg (%.2f lbm)\n', m_ox, m_ox/lbm_to_kg)
fprintf('Fuel CdA: %.2f mm^2\n', fuel_CdA_int*1e6)

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
plot((time(test_ind)), pom(test_ind)/psi_to_Pa,...
    (time(test_ind)), pft(test_ind)/psi_to_Pa,...
    (time(test_ind)), pcc(test_ind)/psi_to_Pa)
xlabel('Time (s)');
ylabel('Pressures (psi)');
legend({'Ox. Manifold','Fuel Tank','Combustion Chamber'})

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
% ExportData(time(test_ind),pft(test_ind),pom(test_ind),pot(test_ind),...
%     zeros(size(time(test_ind))),ft(test_ind),pcc(test_ind))

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

function smoothed = SimpleFilter(time,values,tau,type)
dt = median(diff(time));
Wn = dt/(tau);
values(isnan(values)) = 0;
[b, a] = butter(1,Wn,type);
smoothed = filter(b, a, values);
shift_ind = round(1/(4*Wn));
smoothed(1:(end-shift_ind)) = smoothed((1+shift_ind):end);
end

function smoothed = DiffFilter(time,values,tau)
dt = median(diff(time));
dn = tau/dt;
smoothed = (circshift(values,round(0.5*dn)) - circshift(values,-round(0.5*dn)))./...
     (circshift(time,round(0.5*dn)) - circshift(time,-round(0.5*dn)));
end