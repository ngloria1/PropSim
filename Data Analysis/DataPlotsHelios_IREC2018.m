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
content = fopen( '..\Test Data\IREC_2018.txt', 'r' ) ;
data1 = textscan(content, '%s %s %f %f %f %f %f %f %f', 'Delimiter',', \t', 'HeaderLines', 1 );

%% User Input
t_shutoff_fuel = 9.95; % s, time at shutoff of fuel
rho_fuel = 795; % kg/m^3
gamma_fuelpress = 1.4; % N2
V_ullage_fuel = 0.88e-3; % m^3
pft_start = 725*psi_to_Pa; % psi, after valve open
pft_shutoff = 285*psi_to_Pa; % psi
A_star = pi/4*(2.388e-2)^2; % m^2
m_ox = 4.95; % oxidizer mass expended, kg
dm_fuel = 1e-3*rho_fuel; % fuel mass expended, kg

c_star_theo = 1569; % m/s
Isp_theo = 2079; % m/s
C_f_theo = 1.326; 

%% Import data and convert to SI units
test_time1 = FixTimeOverflow(data1{1,3} / 10^6);
pos = data1{1,4}*psi_to_Pa;
pft = data1{1,5}*psi_to_Pa;
pot = data1{1,6}*psi_to_Pa;

%% Generate Plots
% Test Time Plots
figure()
plot(test_time1, pos/psi_to_Pa,...
     test_time1, pft/psi_to_Pa,...
     test_time1, pot/psi_to_Pa)
xlabel('Time (s)');
ylabel('Pressures (psi)');
legend({'Ox. Supply Tank','Fuel Tank','Ox. Rocket Tank'})

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