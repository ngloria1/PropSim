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

data_Fuel_1 = readtable("../Test Data/2019_05_21_A.csv");
data_Fuel_2 = readtable("../Test Data/2019_05_21_B.csv");
data = readtable("../Test Data/2019_05_21_D.csv");

pft_fill = data_Fuel_1.x4_FuelTank*psi_to_Pa; % Pressure Fuel Tank - fill process
pft_flow = data_Fuel_2.x4_FuelTank*psi_to_Pa; % Pressure Fuel Tank - flowing
pot = data.x4_OxTank*psi_to_Pa; % Pressure Ox Tank
pos = data.x1_NitrousSupply*psi_to_Pa; % Pressure Ox Supply


%% nitrous
Fs = 200; %samples second

time_fuel_fill = linspace(0,length(pft_flow)/Fs,length(pft_fill));
time_fuel_flow = linspace(0,length(pft_flow)/Fs,length(pft_flow));
time_ox = linspace(0,length(pot)/Fs,length(pot));

figure(1)
plot(time_fuel_fill,pft_fill/psi_to_Pa);
xlabel('Time [s]')
ylabel('Fuel Tank Pressure [psi]')
title('Fuel Fill Process')

figure(2)
plot(time_fuel_flow,pft_flow/psi_to_Pa);
xlabel('Time [s]')
ylabel('Fuel Tank Pressure [psi]')
title('Fuel Flow')

figure(3)
plot(time_ox,pot/psi_to_Pa,time_ox,pos/psi_to_Pa);
xlabel('Time [s]')
ylabel('Pressure [psi]')
title('Oxidizer Fill Process')
legend({'Rocket Oxidizer Tank','Supply Oxidizer Pressure'})

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