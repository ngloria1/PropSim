classdef Gas
    %GAS Instance of a gas
    %   Includes properties of a gas
    
    properties
        % molecular mass [kg/mol]
        molecular_mass
        % specific heats [J/kg*K]
        c_v
    end
    
    methods
        function constant = R_specific(obj)
            R_u = 8.3144598; % J/K*mol
            constant = R_u/obj.molecular_mass;
        end
        function specific_heat_pressure = c_p(obj)
            specific_heat_pressure = obj.c_v + obj.R_specific;
        end
        function gamma = gamma(obj)
            gamma = obj.c_p/obj.c_v;
        end
    end
    
end

