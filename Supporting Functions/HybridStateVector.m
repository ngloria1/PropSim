classdef HybridStateVector < StateVector
    %HybridStateVector Object representing the state of a hybrid rocket
    
    properties
        
        m_fuel % fuel mass, kg
    end
    
    methods
        
        function obj = HybridStateVector(inputs, mode)
            if nargin == 2
                obj.inputs = inputs;
                obj.mode = mode;
            end
        end
        
        function V_cc = V_cc(obj)
            % Calculate combustion chamber volume
            
            V_cc_empty = pi/4*(obj.inputs.d_cc)^2*obj.inputs.length_cc;
            V_f = obj.m_fuel/obj.inputs.fuel.rho;
            V_cc = V_cc_empty - V_f;
        end
        
        function rad_port = rad_port(obj)
            % Calculate radius of grain port
            
            V_fuel = obj.m_fuel/obj.inputs.fuel.rho;
            rad_port = sqrt(1/4*obj.inputs.fuel.grain_od^2 - ...
                V_fuel/(pi*obj.inputs.fuel.grain_length));
        end
        
        function column_vector = ColumnVector(obj)
            % Create output column vector
            column_vector = [obj.m_lox;
                obj.m_gox;
                obj.m_oxtank_press;
                obj.m_oxpresstank;
                obj.T_oxtank;
                obj.m_fuel;
                obj.m_cc;
                obj.M_cc;
                obj.gamma_cc;
                obj.T_cc;];
            
        end
    end
    
    
    methods(Static)
        function obj = FromColumnVector(x, inputs, mode)
            obj = HybridStateVector(inputs, mode);
            obj.m_lox = x(1);
            obj.m_gox = x(2);
            obj.m_oxtank_press = x(3);
            obj.m_oxpresstank = x(4);
            obj.T_oxtank = x(5);
            obj.m_fuel = x(6);
            obj.m_cc = x(7);
            obj.M_cc = x(8);
            obj.gamma_cc = x(9);
            obj.T_cc = x(10);
            obj.N2O_properties =  N2O_Properties(obj.T_oxtank);
        end
    end
    
end