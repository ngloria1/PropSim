classdef LiquidStateVector < StateVector
    %LiquidStateVector Object representing the state of a liquid rocket
    
    properties
        T_fueltank_press % fuel tank temperature
        m_fueltank_press  % mass of pressurant in oxidizer tank
        m_fuelpresstank  % mass of pressurant in oxidizer pressurant tank
        m_fuel % fuel mass, kg
    end
    
    methods
        
        function obj = LiquidStateVector(inputs, mode)
            if nargin == 2
                obj.inputs = inputs;
                obj.mode = mode;
            end
        end
        
        function p_fueltank = p_fueltank(obj)
            % Calculate fuel tank pressure
            p_fueltank = obj.m_fueltank_press*...
                obj.inputs.fuel_pressurant.gas_properties.R_specific*...
                obj.T_fueltank_press/obj.V_fuel_ullage;
        end
        
        function p_fuelpresstank = p_fuelpresstank(obj)
            % Calculate pressure in fuel pressurant tank
            m_press_initial = obj.inputs.fuel_pressurant.storage_initial_pressure*...
                obj.inputs.fuel_pressurant.tank_volume/(...
                obj.inputs.fuel_pressurant.gas_properties.R_specific.*...
                obj.inputs.T_amb);
            p_fuelpresstank = obj.inputs.fuel_pressurant.storage_initial_pressure*...
                (obj.m_fuelpresstank/m_press_initial)^...
                obj.inputs.fuel_pressurant.gas_properties.gamma;
        end
        
        function T_fuelpresstank = T_fuelpresstank(obj)
            % Calculate temperature in oxidizer pressurant tank
            T_fuelpresstank = obj.inputs.T_amb*...
                (obj.p_fuelpresstank/obj.inputs.fuel_pressurant.storage_initial_pressure)^...
                ((obj.inputs.fuel_pressurant.gas_properties.gamma-1)/...
                obj.inputs.fuel_pressurant.gas_properties.gamma);
        end
        
        function V_fuel = V_fuel(obj)
            % Calculate volume of fuel
            V_fuel = obj.m_fuel / obj.inputs.fuel.rho;
        end
        
        function V_fuel_ullage = V_fuel_ullage(obj)
            % Calculate ullage volume in fuel tank
            V_fuel_ullage = obj.inputs.fuel.V_tank - obj.V_fuel;
        end
        
        function V_cc = V_cc(obj)
            % Calculate combustion chamber volume
            V_cc = pi/4*(obj.inputs.d_cc)^2*obj.inputs.length_cc;
        end
        
        function column_vector = ColumnVector(obj)
            % Create output column vector
            column_vector = [obj.m_lox;
                obj.m_gox;
                obj.m_oxtank_press;
                obj.m_oxpresstank;
                obj.m_fueltank_press;
                obj.m_fuelpresstank;
                obj.T_oxtank;
                obj.T_fueltank_press;
                obj.m_fuel;
                obj.m_cc;
                obj.M_cc;
                obj.gamma_cc;
                obj.T_cc;];
        end
    end
    
    
    methods(Static)
        function obj = FromColumnVector(x, inputs, mode)
            obj = LiquidStateVector(inputs, mode);
            obj.m_lox = x(1);
            obj.m_gox = x(2);
            obj.m_oxtank_press = x(3);
            obj.m_oxpresstank = x(4);
            obj.m_fueltank_press = x(5);
            obj.m_fuelpresstank = x(6);
            obj.T_oxtank = x(7);
            obj.T_fueltank_press = x(8);
            obj.m_fuel = x(9);
            obj.m_cc = x(10);
            obj.M_cc = x(11);
            obj.gamma_cc = x(12);
            obj.T_cc = x(13);
            obj.N2O_properties =  N2O_Properties(obj.T_oxtank);
        end
    end
    
end