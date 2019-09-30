function grad = FindMassGradient(constraints, values, delta, results_prev)
    %% Sample base point
    % Create design, goal structures
    [goal, design] = CreateStructs(constraints, values);
    
    base_inputs = DesignLiquid(results_prev, goal, design, false);
    central_value = CalcMass(base_inputs);
    
    % Sample gradient
    grad = zeros(1, length(values));
    for j = 1:length(values)
        step = zeros(size(values));
        step(j) = values(j)*delta;
        
        [goal, design] = CreateStructs(constraints, values + step);
        
        grad(j) = (CalcMass(DesignLiquid(base_inputs, goal, design, false)) - central_value)/delta;
    end
    grad = grad/norm(grad);
end

function [goal, design] = CreateStructs(constraints, values)

    % Create design, goal structures
    goal.max_thrust = constraints(2);
    goal.OF = constraints(3);
    goal.total_impulse = constraints(1);
    goal.min_fuel_dp = constraints(4); % min dp as % of tank pressure
    goal.min_ox_dp = constraints(5); % min dp as % of tank pressure
    goal.ox_to_fuel_time = values(1); % ratio of liquid oxidizer flow time to liquid fuel flow time

    design.p_tanks = values(2);
    design.ox_ullage = constraints(6);
    design.exp_ratio = values(3);
    
end