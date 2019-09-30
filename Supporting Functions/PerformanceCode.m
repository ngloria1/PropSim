function record = PerformanceCode(inputs, mode, test_data, options)
% Take inputs defining motor characteristics, take test data, and calculate
% and plot performance. 
%   INPUTS:
%       - inputs: see Integration.m
%       - mode: see Integration.m
%       - test_data: structure of data related to test data for plotting
%       purposes
%           - test_plots_on: 1 for plot test data, 0 for nothing
%           - test_data_file: filename of test data within "Test Data"
%           directory
%               - All data should be in physical values (not voltages) in
%               SI units of m, kg, K, s, Pa, N, etc.
%           - t_offset: time at which ignition occurs in test data]
%       - options: set of options
%           - t_final: simulation time limit
%           - dt: output time resolution
%           - output_on: true for plots, false for no plots
%   OUTPUTS:
%       - Isp: specific impulse [s]
%       - time: time vector over burn for F_thrust
%       - F_thrust: vector of thrust over burn
%       - F_thrust_RASAERO: vector of thrust points used for RASAero thrust
%       curve output
%       - impulse: total impulse
%       - pressure_drop: average pressure drop over burn duration
%       - F_thrust_RASAERO.txt file in "Outputs" with thrust curve that can 
%       be adap_oxtanked for RASAero or open rocket. 

% To Do:
%   - Change to higher-order integration method (e.g. Runge-Kutta 4th order
%   or matlab ODE45 function
%   - Add feedline fricitional losses to two-phase injector flow calculator
%   - Investigate differences between simulation and test data for
%   supercharged oxidizer tank case

%% Constants
R_u = 8.3144621; % universal gas constant [J/mol*K]
g_0 = 9.80665; % standard gravitational constant [m/s^2]

% Unit Conversion
psi_to_Pa = 6894.75729; % 1 psi in Pa
in_to_m = 0.0254; % 1 in in m
lbf_to_N = 4.44822162; % 1 lbf in N
atm_to_Pa = 101325; % 1 atm in Pa

%% Integration Parameters
default_options.t_final  =  60;    % Integration time limit
default_options.dt      = 0.01;  % Timestep [s]
default_options.output_on = true;
if nargin < 4
    options = default_options;
else
    if ~isfield(options, 't_final')
        options.t_final = default_options.t_final;
    end
    if ~isfield(options, 'dt')
        options.dt = default_options.dt;
    end
    if ~isfield(options, 'output_on')
        options.output_on = default_options.output_on;
    end
end
tspan = 0:options.dt:options.t_final;

%% Calculate initial properties of the nitrous in the tank
N2O = N2O_Properties(inputs.ox.T_tank);

%Our integration variables are oxidizer mass and liquid oxidizer volume
Mox = N2O.rho_l*(inputs.ox.V_l) + N2O.rho_g*(inputs.ox.V_tank - inputs.ox.V_l);
if options.output_on
    fprintf('Initial oxidizer mass: %.2f kg\n', Mox);
end

%% Create Injector Area Function
function [F] = Throttle(t)
    F = zeros(length(t));
    for i = 1:length(t)
        if t(i) > inputs.dt_valve_open
            F(i) = 1;
        else
            F(i) = (t(i)/inputs.dt_valve_open);
        end
    end
end
inputs.Throttle = @Throttle;

%% Run Integration Function
tic
[time, record] = ...
    Integration(inputs,mode,tspan);

F_thrust = record.F_thrust;
p_cc = record.p_cc;
p_oxtank = record.p_oxtank;
p_oxpresstank = record.p_oxpresstank;
p_fueltank = record.p_fueltank;
p_fuelpresstank = record.p_fuelpresstank;
p_oxmanifold = record.p_oxmanifold;
T_oxtank = record.T_oxtank;
T_cc = record.T_cc;
area_core = record.area_core;
OF = record.OF_i;
gamma_ex = record.gamma_ex;
m_dot_ox = record.m_dot_ox;
m_dot_fuel = record.m_dot_fuel;
p_crit = record.p_crit;
m_dot_ox_crit = record.m_dot_ox_crit;
M_e = record.M_e;
p_exit = record.p_exit;
p_shock = record.p_shock;
if options.output_on
    toc
end

%% Plot Results
if options.output_on
    if test_data.test_plots_on
        % Load Imported Data for comparison
        [test_time, pft, pom, pot, we, ft, pcc] = LoadDataVars(test_data.test_data_file, test_data.t_offset);
    else
        test_time = [];
        pft = [];
        pom = [];
        pot = [];
        we = [];
        ft = [];
        pcc = [];
    end

    h_figure_1 = figure;
    h_tabgroup = uitabgroup(h_figure_1);
    tab1 = uitab(h_tabgroup,'Title','Thrust');
    tab2 = uitab(h_tabgroup,'Title','Pressures');
    tab3 = uitab(h_tabgroup,'Title','Pressure Drop');
    tab4 = uitab(h_tabgroup,'Title','Temperatures');
    tab5 = uitab(h_tabgroup,'Title','Fuel');
    tab6 = uitab(h_tabgroup,'Title','Oxidizer Mass Flux');
    tab7 = uitab(h_tabgroup,'Title','Performance');
    tab8 = uitab(h_tabgroup,'Title','Nozzle');

    if mode.combustion_on
        axes('parent',tab1)
        plot(time, F_thrust./lbf_to_N,...
        test_time,ft./lbf_to_N)
        xlabel('Time (s)')
        ylabel('Thrust Force (lbf)')
    end

    axes('parent',tab2)
    if mode.combustion_on
        subplot(2,1,1)
        plot(time,p_cc(1:length(time))./psi_to_Pa,...
        test_time,pcc./psi_to_Pa)
        xlabel('Time (s)')
        ylabel('Chamber Pressure (psi)')
        if isempty(test_time)
            legend('Simulation')
        else
            legend('Simulation', 'Measured')
        end
    end

    subplot(2,1,2)
    legend_str = {};
    if inputs.ox_pressurant.active && ...
        inputs.ox_pressurant.tank_volume > 0
        plot(time, p_oxpresstank(1:length(time))./psi_to_Pa)
        hold on
        legend_str(end+1) = {'Oxidizer Pressurant Tank (sim)'};
    end
    if strcmp(mode.type, 'liquid') && inputs.fuel_pressurant.active && ...
        inputs.fuel_pressurant.tank_volume > 0
        plot(time, p_fuelpresstank(1:length(time))./psi_to_Pa)
        hold on
        legend_str(end+1) = {'Fuel Pressurant Tank (sim)'};
    end
    plot(time,p_oxtank(1:length(time))./psi_to_Pa, ...
        time, p_oxmanifold(1:length(time))./psi_to_Pa,...
        time,p_fueltank(1:length(time))./psi_to_Pa,...
        test_time,pot./psi_to_Pa,...
        test_time,pom./psi_to_Pa,...
        test_time,pft./psi_to_Pa)
    xlabel('Time (s)')
    ylabel('Tank Pressure (psi)')
    legend_str(end+1:end+3) = {'Oxidizer Tank (sim)', 'Oxidizer Manifold (sim)'...
        'Fuel Tank (sim)'};
    if ~isempty(test_time)
        legend_str(end+1:end+3) = {'POT (test)','POM (test)','PFT (test)'};
    end
    legend(legend_str)

    axes('parent',tab3)
    legend_str = {'Oxidizer (sim)', 'Oxidizer Critical (sim)'};
    if strcmp(mode.type, 'liquid')
        n_plots = 2;
    else
        n_plots = 1;
    end
    subplot(n_plots,1,1)
    plot(time, 100*record.ox_pressure_drop, time, 100*(1-p_crit./p_oxtank))
    if ~isempty(test_time)
        legend_str(end+1) = {'Measured'};
    end
    xlabel('Time (s)')
    ylabel({'Pressure Drop across Injector', '[% of Tank Pressure]'})
    legend(legend_str)

    if strcmp(mode.type, 'liquid')
        subplot(n_plots,1,2)
        hold on
        plot(time, 100*record.fuel_pressure_drop)
        xlabel('Time (s)')
        ylabel({'Pressure Drop across Injector', '[% of Tank Pressure]'})
        legend('Fuel (sim)')
    end

    axes('parent',tab4)
    if mode.combustion_on
        subplot(2,1,1)
        plot(time,T_cc(1:length(time)))
        xlabel('Time (s)')
        ylabel('Chamber Temperaure [k]')
        legend('Simulation')
    end

    subplot(2,1,2)
    plot(time,T_oxtank(1:length(time))-273.15)
    xlabel('Time (s)')
    ylabel('Tank Temperature (C)')
    if isempty(test_time)
        legend('Simulation')
    else
        legend('Simulation', 'Measured (Top)', 'Measured (Bottom)')
    end

    axes('parent',tab5)
    subplot(2,2,1)
    plot(time,sqrt(area_core/pi)*2/in_to_m)
    xlabel('Time (s)')
    ylabel('Grain Port Diameter (in)')
    legend('Simulation')

    subplot(2,2,2)
    plot(time,m_dot_ox,time,m_dot_ox_crit)
    xlabel('Time (s)')
    ylabel('Oxidizer Mass Flow Rate [kg/s]')
    legend({'Simulation','Simulation Choked'})

    subplot(2,2,3)
    plot(time,m_dot_fuel)
    xlabel('Time (s)')
    ylabel('Fuel Mass Flow Rate [kg/s]')
    legend('Simulation')

    subplot(2,2,4)
    plot(time,OF)
    xlabel('Time (s)')
    ylabel('OF Ratio')
    legend('Simulation')

    if mode.combustion_on
        axes('parent',tab6)
        plot(time, m_dot_ox./area_core)
        xlabel('Time (s)')
        ylabel('Grain Core Oxidizer Mass Flux (kg/m^{2}s)')
        legend('Simulation')
    end

    axes('parent',tab7)
    if mode.combustion_on
        subplot(2,1,1)
        plot(time, record.Isp_i, time, record.c_star_i)
        xlabel('Time (s)')
        ylabel('Velocity [m/s]')
        legend('I_{sp}', 'C*')

        subplot(2,1,2)
        plot(time, record.c_f_i)
        xlabel('Time (s)')
        ylabel('Coefficient of Thrust (C_f) []')
    end

    axes('parent',tab8)
    if mode.combustion_on
        subplot(2,1,1)
        plot(time, p_exit*atm_to_Pa^-1, time, p_shock*atm_to_Pa^-1)
        xlabel('Time (s)')
        ylabel('Nozzle Pressure (atm)')
        legend('Sim. Exit Pressure', 'Sim. Back Pressure Threshold for Normal Shock')

        subplot(2,1,2)
        plot(time, M_e)
        xlabel('Time (s)')
        ylabel('Exit Mach Number ()')
        legend('Simulation')
    end

    %Use trapezoidal integration
    impulse = trapz(time, F_thrust);
    Mox_initial = record.m_ox(1);
    Mox = record.m_ox(1) - record.m_ox(end);

    Mfuel_initial = record.m_fuel(1);
    Mfuel = record.m_fuel(1) - record.m_fuel(end);

    fprintf('Pressurant Mass: %.3f kg\n', record.m_press)
    fprintf(['Impulse: %.2f kN*s\t\t\nOxidizer Mass Spent: '...
        '%.2f kg\t\tOxidizer Mass Remaining: %.2f kg\nFuel Mass Spent: %.2f kg\t\t' ...
        'Fuel Mass Remaining: %.2f kg\n OF ratio: %.2f \n']...
        , impulse/1000, Mox, Mox_initial-Mox, Mfuel, Mfuel_initial - Mfuel, ...
        Mox/Mfuel);
    fprintf('Isp: %.1f s\t\tC*: %.0f m/s\t\tC_f: %.2f\n', ...
        record.Isp/g_0, record.c_star, record.c_f)

    % F_thrust_RASAERO must be less than or equal to 32 entries
    % Inputs
    num_entries = 30;
    motor_designation = 'I';           
    length_motorcase = 2438.4; %casing length in millimeters
    delay = 'P'; %always P for plugged (no ejection charge)
    dry_weight = 11.94; %kg
    diameter = 143; %motor casing diameter in millimeters
    manufacturer = 'SSI';
    switch impulse
        case (impulse < 10200)
            motor_class = 'M';
        case (impulse >= 10200) && (impulse < 20500)
            motor_class = 'N';
        otherwise
            motor_class = 'O';
    end

    engine_name = sprintf('%s%.0f-%s', motor_class, 100 * round(impulse/100), motor_designation);
    propellant_weight = Mox + Mfuel; %kg
    tot_weight = dry_weight + propellant_weight; %kg

    F_thrust_RASAERO = zeros(num_entries,2);
    F_thrust_RASAERO(1:(num_entries),1) = time(round(linspace(2,length(F_thrust),num_entries)))';
    F_thrust_RASAERO(1:(num_entries),2) = F_thrust(round(linspace(2,length(F_thrust),num_entries)))';
    F_thrust_RASAERO(num_entries,:) = [time(end), 0];

    fid = fopen('./Outputs/F_thrust_RASAERO.txt','w');
    fprintf(fid, '; Name diameter(mm) Length(mm) delay propellant_weight(kg) mass(kg)\n');
    fprintf(fid, '%s %.0d %.0f %s %.2f %.2f %s\n', engine_name, diameter, length_motorcase, delay, propellant_weight, tot_weight, manufacturer);
    fprintf(fid, '%.3f %.3f\n', F_thrust_RASAERO');
    fclose(fid);
end
end

function [test_time, pft, pom, pot, we, ft, pcc] = LoadDataVars(filename, t_offset)
    test_time = 0;
    pft = 0;
    pom = 0;
    pot = 0;
    we = 0;
    ft = 0;
    pcc = 0;
    load(['.\Test Data\' filename])
    test_time = test_time + t_offset;
    if length(test_time) == 1
        error('Variable ''test_time'' not present in imported data.\n');
    end
    
%     Create emp_oxtanky vectors for missing data
%     n_indices = length(test_time);
%     if length(force_thrust) == 1
%         fprintf('Warning: ''force_thrust'' missing!\n');
%         force_thrust = zeros(1,n_indices);
%     end
%     if length(tank_top_pressure) == 1
%         fprintf('Warning: ''tank_top_pressure'' missing!\n');
%         tank_top_pressure = zeros(1,n_indices);
%     end
%     if length(tank_bottom_pressure) == 1
%         fprintf('Warning: ''tank_bottom_pressure'' missing!\n');
%         tank_bottom_pressure = zeros(1,n_indices);
%     end
%     if length(tank_top_temperature) == 1
%         fprintf('Warning: ''tank_top_temperature'' missing!\n');
%         tank_top_temperature = zeros(1,n_indices);
%     end
%     if length(tank_bottom_temperature) == 1
%         fprintf('Warning: ''tank_bottom_temperature'' missing!\n');
%         tank_bottom_temperature = zeros(1,n_indices);
%     end
%     if length(combustion_chamber_pressure) == 1
%         fprintf('Warning: ''combustion_chamber_pressure'' missing!\n');
%         combustion_chamber_pressure = zeros(1,n_indices);
%     end
%     if length(pressure_drop) == 1
%         fprintf('Warning: ''pressure_drop'' missing!\n');
%         pressure_drop = zeros(1,n_indices);
%     end
end