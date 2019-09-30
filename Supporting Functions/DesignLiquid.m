function [inputs] = DesignLiquid(initial_inputs, goal, design, output_on)
%%DesignLiquid Generates a liquid engine design matching the performance
%%specifications and constraints provided.
%   Inputs:
%       - goal.max_thrust: initial thrust required
%       - goal.OF: OF ratio desired
%       - goal.total_impulse: total impulse desired
%       - goal.min_fuel_dp: minimum pressure drop from fuel tank to
%       injector as proportion of fuel tank pressure
%       - goal.min_ox_dp: minimum pressure drop from oxidizer tank to
%       injector as proportion of oxidizer tank pressure
%       - goal.ox_to_fuel_time: ratio of liquid oxidizer flow time to liquid fuel flow time
%       - design.p_tanks: initial tank pressures
%       - design.ox_ullage: initial oxidizer tank volume ullage fraction
%       - design.exp_ratio: design expansion ratio
%   Outputs:
%       - inputs: engine design parameters (see PerformanceCode.m)

format long

%% Get Parameter Field Names
parameter_field_names = fieldnames(goal);
parameter_labels = {'Thrust','OF','Total Impulse','Min. Fuel Pressure Drop',...
    'Min. Ox. Pressure Drop','Ox. to Fuel Drain Time Ratio'};

%% Options
test_data.test_plots_on = 0; % Import tests data and plot against simulation data
test_data.test_data_file = ''; % File from which to import test data
test_data.t_offset = 0; % Time offset of test data wrt simulation data [s]

% 1: simulate combustion (hot fire), 0: no combustion (cold flow)
mode.combustion_on = 1;
% 1: simulate flight conditions (i.e. acceleration head), 0: ground test 
% conditions
mode.flight_on = 0;
% 'hybrid' for solid fuel, liquid oxidizer, 'liquid' for liquid fuel and
% oxidizer
mode.type = 'liquid';

%% Run Performance Code
max_iter = 100;
param_tol = 0.005;
% Enforce design parameters
inputs = initial_inputs;
inputs.ox.T_tank = N2O_FindT(design.p_tanks);
inputs.ox.V_tank = inputs.ox.V_l/(1-design.ox_ullage);
inputs.fuel_pressurant.set_pressure = design.p_tanks;
inputs.exp_ratio = design.exp_ratio;

% Initialize residual plot
if output_on
    figure
    resid_axes = axes;
    figure
    value_axes = axes;
end
errors = cell2struct(cell(length(parameter_field_names),1),parameter_field_names);
for ii = 1:length(parameter_field_names)
    errors.(parameter_field_names{ii}) = NaN*ones(max_iter,1);
end
parameters = errors;

for ii = 1:max_iter
    options.output_on = true;
    options.dt = 0.005; % small dt for increased accuracy
    record = PerformanceCode(inputs, mode, test_data, options);
    dt_thrust_filter = 0.1;
    [max_thrust,m_dot_max,p_cc_max] = ...
        FindMaxThrust(record,dt_thrust_filter);
    
    [t_fuel_liq, t_ox_liq] = FindTLiqs(record);
    
    [parameters_ii,errors_ii] = CalculateParamErrors(goal, max_thrust, ...
        t_ox_liq/t_fuel_liq, record, parameter_field_names);
    for jj = 1:length(parameter_field_names)
        errors.(parameter_field_names{jj})(ii) = ...
            errors_ii.(parameter_field_names{jj});
        parameters.(parameter_field_names{jj})(ii) = ...
            parameters_ii.(parameter_field_names{jj});
    end
    
    if output_on
        % Plot residuals
        cla(resid_axes)
        for jj = 1:length(parameter_field_names)
            semilogy(resid_axes,1:ii,errors.(parameter_field_names{jj})(1:ii))
            hold(resid_axes,'on')
        end
        title(resid_axes,'Design Parameter Errors')
        xlabel(resid_axes,'Iteration #')
        ylabel(resid_axes,'Relative Error')
        legend(resid_axes,parameter_labels)

        cla(value_axes)
        for jj = 1:length(parameter_field_names)
            plot(value_axes,1:ii,parameters.(parameter_field_names{jj})(1:ii)...
                /goal.(parameter_field_names{jj}))
            hold(value_axes,'on')
        end
        title(value_axes,'Design Parameters')
        xlabel(value_axes,'Iteration #')
        ylabel(value_axes,'Relative Values')
        legend(value_axes,parameter_labels)
        drawnow
    end
    
    % Terminate iteration once desired results achieved
    parameter_converged = zeros(length(parameter_field_names),1);
    % Check that value is below tolerance
    for jj = 1:length(parameter_field_names)
            parameter_converged(jj) = errors.(parameter_field_names{jj})(ii) < param_tol;
    end
    if all(parameter_converged)
        break
    end
    
    % Ox/fuel injector area distribution
    ox_to_fuel_time_change = goal.ox_to_fuel_time/(t_ox_liq/t_fuel_liq);
    % Keep mass flow rate constant
    OF_time_fuel_change = 1 + (ox_to_fuel_time_change-1)*(goal.OF/(1+goal.OF));
    OF_time_ox_change = 1 - (ox_to_fuel_time_change-1)*(1/(1+goal.OF));
    
    % Tank sizing
    impulse_change = goal.total_impulse/record.impulse;
    % Disable OF changes unless OF time is close to converged. OF changes
    % assumes that all fuel/oxidizer is being spent. This assumption is not
    % valid if the engine cuts out with excess fuel or oxidizer. 
    if errors.ox_to_fuel_time(ii) < 0.01
        OF_change = goal.OF/record.OF;
    else
        OF_change = 1;
    end
    % Keep mass flow rate constant, while shifting OF
    OF_fuel_change = 1 - (OF_change-1)*(goal.OF/(1+goal.OF));
    OF_ox_change = 1 + (OF_change-1)*(1/(1+goal.OF));
    inputs.ox.V_tank = inputs.ox.V_tank*impulse_change*OF_ox_change;
    inputs.ox.V_l = inputs.ox.V_l*impulse_change*OF_ox_change;
    inputs.fuel.V_tank = inputs.fuel.V_tank*impulse_change*OF_fuel_change;
    inputs.fuel.V_l = inputs.fuel.V_l*impulse_change*OF_fuel_change;
    
    % Thrust sizing
    ox_pressure_drop_change = goal.min_ox_dp/parameters.min_ox_dp(ii);
    fuel_pressure_drop_change = goal.min_fuel_dp/parameters.min_fuel_dp(ii);
    % Account for oxidizer pressure drop change by raising combustion
    % chamber pressure
    p_cc_change = (1-goal.min_ox_dp)/(1-parameters.min_ox_dp(ii));
    
    V_fuel_change = fuel_pressure_drop_change^(2/3); 
        % 2/3 constant empirical fit
    m_dot_change = goal.max_thrust/max_thrust;
    
    A_star_change = 1/p_cc_change*m_dot_change;
    
    inputs.ox.injector_area = inputs.ox.injector_area*m_dot_change...
        *OF_time_ox_change*OF_ox_change;
    inputs.fuel.injector_area = inputs.fuel.injector_area*m_dot_change...
        *OF_time_fuel_change*OF_fuel_change;
    
    inputs.d_throat = A_star_change^0.5*inputs.d_throat;
    
    inputs.fuel.V_tank = max(inputs.fuel.V_tank*V_fuel_change,inputs.fuel.V_l*1.01);
end

if output_on
    PrintResults(inputs);
    PerformanceCode(inputs, mode, test_data);
end

end

function [max_thrust,m_dot_max,p_cc_max] = ...
    FindMaxThrust(record,dt_filter)
%%FindMaxThrust Find the maximum thrust on the thrust curve, filtering out
%%spikes of duration less than the filter time constant.

dn_thrust_filter = ceil(dt_filter/mean(diff(record.time)));
a = 1;
b = 1/dn_thrust_filter*ones(dn_thrust_filter,1);
filtered_thrust = filter(b,a,record.F_thrust);
filtered_m_dot = filter(b,a,record.m_dot_prop);
filtered_p_cc = filter(b,a,record.p_cc);
[~,max_ind] = max(filtered_thrust);    
max_thrust = filtered_thrust(max_ind);
m_dot_max = filtered_m_dot(max_ind);
p_cc_max = filtered_p_cc(max_ind);
end

function [t_fuel_liq, t_ox_liq] = FindTLiqs(record)
%%FindTLiqs Find liquid drain times, given performance record.

if ~isnan(record.t_fuel_liq)
    t_fuel_liq = record.t_fuel_liq;
else
    t_fuel_liq = interp1(record.m_fuel,record.time,0,'linear','extrap');
end
if ~isnan(record.t_ox_liq)
    t_ox_liq = record.t_ox_liq;
else
    dm_tol = 1e-3;
    if record.m_lox(end) < dm_tol
        t_ox_liq = record.time(find(record.m_lox < dm_tol,1));
    else
        t_ox_liq = interp1(record.m_lox,record.time,0,'linear','extrap');
    end
end
end

function [parameters,errors] = CalculateParamErrors(goal, max_thrust, ...
    ox_to_fuel_time, record, field_names)
%%CalculateParamErrors Calculate the relative errors in goal parameters

parameters.max_thrust = max_thrust;
parameters.total_impulse = record.impulse;
parameters.OF = record.OF;
parameters.ox_to_fuel_time = ox_to_fuel_time;
parameters.min_ox_dp = min(record.ox_pressure_drop);
parameters.min_fuel_dp = min(record.fuel_pressure_drop);

for ii = 1:length(field_names)
    errors.(field_names{ii}) = ...
        abs(parameters.(field_names{ii}) - goal.(field_names{ii}))/...
        goal.(field_names{ii});
end
end

function PrintResults(inputs)
%%PrintResults Print design characteristics as a result of iteration.


% Unit Conversion
psi_to_Pa = 6894.75729; % 1 psi in Pa
in_to_m = 0.0254; % 1 in in m
mm_to_m = 1e-3; % 1 mm in m
lbf_to_N = 4.44822162; % 1 lbf in N
lbm_to_kg = 0.453592; % 1 lbm in kg
atm_to_Pa = 101325; % 1 atm in Pa
L_to_m3 = 1e-3; % 1 L in m^3

fprintf('Converged Results:\n')
fprintf('Nozzle Throat Diameter: %.3f cm\n', inputs.d_throat*1e2)
fprintf('Fuel Injector CdA: %.3g mm^2\n', inputs.fuel.injector_area*(1e3)^2)
fprintf('Oxidizer Injector CdA: %.3g mm^2\n', inputs.ox.injector_area*(1e3)^2)
fprintf('Fuel Tank Size: %.3g L tank, %.3g L liquid\n', ...
    inputs.fuel.V_tank/L_to_m3, inputs.fuel.V_l/L_to_m3)
fprintf('Oxidizer Tank Size: %.3g L tank, %.3g L liquid\n', ...
    inputs.ox.V_tank/L_to_m3, inputs.ox.V_l/L_to_m3)
fprintf('Fuel Tank Initial Pressure: %.3g psi\n', ...
    inputs.fuel_pressurant.set_pressure/psi_to_Pa)
fprintf('Oxidizer Tank Initial Temperature: %.3g K\n', ...
    inputs.ox.T_tank)
end