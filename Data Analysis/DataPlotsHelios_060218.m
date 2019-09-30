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

%content = fopen( 'C:\Users\James\Desktop\AA284B\POM-PFT-POT-Weight-Thrust.txt', 'r' ) ;
content = fopen( '..\Test Data\6_02_18_pressures.txt', 'r' ) ;
data1 = textscan(content, '%s %s %f %f %f %f %f %f %f', 'Delimiter',', \t', 'HeaderLines', 1 );
content2 = fopen( '..\Test Data\6_02_18_cc.txt', 'r' ) ;
data2 = textscan(content2, '%s %s %f %f %f %f', 'Delimiter',', \t', 'HeaderLines', 0 );

%% User Input
t_shutoff_fuel = 9.95; % s, time at shutoff of fuel
rho_fuel = 795; % kg/m^3
gamma_fuelpress = 1.4; % N2
V_ullage_fuel = 0.88e-3; % m^3
pft_start = 750*psi_to_Pa; % psi, after valve open
pft_shutoff = 350*psi_to_Pa; % psi
pft_empty = 278*psi_to_Pa; % psi
A_star = pi/4*(2.388e-2)^2; % m^2
m_ox = 4.95; % oxidizer mass expended, kg
% dm_fuel = 1e-3*rho_fuel; % fuel mass expended, kg
dm_fuel = rho_fuel * V_ullage_fuel*((pft_start./pft_shutoff).^(1/gamma_fuelpress)-1);
m_fuel_loaded = rho_fuel * V_ullage_fuel*((pft_start./pft_empty).^(1/gamma_fuelpress)-1);
dm_fuel = dm_fuel/m_fuel_loaded*(1e-3*rho_fuel);

c_star_theo = 1579; % m/s
Isp_theo = 2079; % m/s
C_f_theo = 1.326; 

%% Import data and convert to SI units
test_time1 = FixTimeOverflow(data1{1,3} / 10^6);
pom = data1{1,4}*psi_to_Pa;
pft = data1{1,5}*psi_to_Pa;
pot = data1{1,6}*psi_to_Pa;
we = data1{1,7}*lbf_to_N;
ft = data1{1, 8}*lbf_to_N;
test_time2 = FixTimeOverflow(data2{1,3}/ 10^6) + 7.8;
pcc = data2{1,4}*psi_to_Pa;

%% Find burn end and start times
pom_threshold = 50*psi_to_Pa;
pcc_threshold = 50*psi_to_Pa;
tau_filter = 0.1; % s
pom_filt = SimpleFilter(test_time1,pom,tau_filter,'low');
pcc_filt = SimpleFilter(test_time2,pcc,tau_filter,'low');
burn_start_time = test_time1(find(pom_filt>pom_threshold,1)) - 0.05;
burn_end_time = test_time1(find(pom_filt<pom_threshold & test_time1>burn_start_time + 0.1,1)) + 0.05;
burn_start_time2 = test_time2(find(pcc_filt>pcc_threshold,1)) - 0.05;
dt_shift_2 = burn_start_time2 - burn_start_time;

test_time2 = test_time2 - dt_shift_2;

%% Change to common time index
time = sort(unique([test_time1; test_time2]));
pom = interp1(test_time1,pom,time);
pft = interp1(test_time1,pft,time);
pot = interp1(test_time1,pot,time);
we = interp1(test_time1,we,time);
ft = interp1(test_time1,ft,time);
pcc = interp1(test_time2,pcc,time);

dp_fuel = pft - pcc;
dp_ox = pom - pcc;

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

% Tare weight to zero before test start
we_offset = mean(we((- 5 < time) & (time < 0)));
we = we - we_offset;
ft_offset = mean(ft((- 5 < time) & (time < 0)));
ft = ft - ft_offset;

% Take thrust as combination of thrust load cell and removed weight
ft_raw = ft;
ft = ft + we;

%% Filtering
tau_filter = 1.0; % s
pft_filt = SimpleFilter(time,pft,tau_filter,'low');
dp_fuel_filt = SimpleFilter(time,dp_fuel,tau_filter,'low');
pcc_diff = SimpleFilter(time,dp_fuel,tau_filter,'high');

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
c_star_int = (pcc_int*A_star)/(m_ox + dm_fuel_test);
Isp = Impulse/(m_ox + dm_fuel_test);
OF = m_ox/dm_fuel_test;
burn_time = 7.0; %max(time(burn_ind));
ft_avg = Impulse/burn_time;
pcc_avg = pcc_int/burn_time;

fprintf('Calculations:\n')
fprintf('Impulse: %.3g kN*s (%.3g lbf*s)\n', Impulse/1e3, Impulse/lbf_to_N)
fprintf('Burn Time: %.3g s\n', burn_time)
fprintf('Avg. Thrust: %.3g kN (%.3g lbf)\n', ft_avg/1e3, ft_avg/lbf_to_N)
fprintf('Avg. Pcc: %.3g MPa (%.3g psi)\n', pcc_avg/1e6, pcc_avg/psi_to_Pa)
fprintf('Avg. Pcc RMS: %.3g %%\n', pcc_rms_int/pcc_int*1e2)
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
plot((time(test_ind)), dp_fuel(test_ind)./pft(test_ind),...
    (time(test_ind)), dp_ox(test_ind)./pot(test_ind))
xlabel('Time (s)');
ylabel('Pressure Drop (% of tank pressure)');
legend({'Fuel','Oxidizer'})
figure()
plot((time(test_ind)), ft_raw(test_ind)/lbf_to_N,...
    (time(test_ind)), we(test_ind)/lbf_to_N)
xlabel('Time (s)');
ylabel('Force (lbf)');
legend({'Thrust Load Cell','Weight Load Cell'})
figure()
plot((time(test_ind)), ft(test_ind)/lbf_to_N)
xlabel('Time (s)');
ylabel('Thrust (lbf)');

% Performance Parameters
figure()
plot(time(test_ind), C_f(test_ind))
xlabel('Time (s)');
ylabel('Thrust Coefficient ()');

% CdA Plots
figure()
plot(time(time > -1 & time < t_shutoff_fuel + 10),pft_filt(time > -1 & time < t_shutoff_fuel + 10)/psi_to_Pa,...
    time(time > -1 & time < t_shutoff_fuel + 10),pft(time > -1 & time < t_shutoff_fuel + 10)/psi_to_Pa)
xlabel('Time (s)');
ylabel('Pressure Fuel Tank - filtered (psi)');
figure()
plot(time(time > -1 & time < t_shutoff_fuel + 10),m_fuel(time > -1 & time < t_shutoff_fuel + 10))
xlabel('Time (s)');
ylabel('Fuel Mass Consumed (kg)');
figure()
plot(time(time > -1 & time < t_shutoff_fuel + 10),m_dot_fuel(time > -1 & time < t_shutoff_fuel + 10))
xlabel('Time (s)');
ylabel('Fuel Mass Flow Rate (kg/s)');
figure()
plot(time(nominal_flow_ind),fuel_CdA(nominal_flow_ind)*1e6)
ylim([0 2*median(fuel_CdA(nominal_flow_ind))*1e6])
xlabel('Time (s)');
ylabel('Fuel CdA (mm^2)');

%% Export Data
ExportData(time(test_ind),pft(test_ind),pom(test_ind),pot(test_ind),...
    we(test_ind),ft(test_ind),pcc(test_ind))

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