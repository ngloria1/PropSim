classdef StateVector
    %STATEVECTOR Value of engine state, allowing access to variables by
    %name
    
    properties
        inputs
        mode
        
        % Tanks
        m_lox
        m_gox
        m_oxtank_press  % mass of pressurant in oxidizer tank
        m_oxpresstank  % mass of pressurant in oxidizer pressurant tank
        T_oxtank
        N2O_properties
        
        % Combustion chamber
        m_cc  % mass of gas
        M_cc  % molecular mass of gas
        gamma_cc  % ratio of specific heats of gas
        T_cc  % temperature of gas
    end
   
    methods
        
        function obj = StateVector(inputs, mode)
            if nargin == 2
                obj.inputs = inputs;
                obj.mode = mode;
            end
        end
        
        function V_lox = V_lox(obj)
            % Calculate volume of liquid oxidizer
            V_lox = obj.m_lox / obj.N2O_properties.rho_l;
        end
        
        function p_oxtank = p_oxtank(obj)
            % Calculate volume of liquid oxidizer
            p_oxtank = obj.p_gox + obj.p_oxtank_press;
        end
        
        function m_cv = oxtank_m_cv(obj)
            % Thermal capacity of oxidizer tank (constant volume)
            m_cv = (obj.m_lox*obj.N2O_properties.cv_l + ...
                obj.m_gox*obj.N2O_properties.cv_g + ...
                obj.m_oxtank_press*obj.inputs.ox_pressurant.gas_properties.c_v);
        end
        
        function m_cp = oxtank_m_cp(obj)
            % Thermal capacity of oxidizer tank (constant pressure)
            m_cp = (obj.m_lox*obj.N2O_properties.cp_l + ...
                obj.m_gox*obj.N2O_properties.cp_g + ...
                obj.m_oxtank_press*obj.inputs.ox_pressurant.gas_properties.c_p);
        end
        
        function p_gox = p_gox(obj)
            % Calculate oxidizer tank partial pressure of oxidizer
            R_u = 8.3144621; % Universal gas constant [J/mol*K]
            M_n2o = 0.044013; % Molecular mass of nitrous oxide [kg/mol]
            R_n2o = R_u/M_n2o; %  Specific gas constant of nitrous oxide [J/kg*K]
            a_n2o = 0.38828/M_n2o^2; % van der Waal's constant a for N2O [[Pa*(kg/m^3)^-2]]
            b_n2o = 44.15/M_n2o*10^-6; % van der Waal's constant a for N2O [m^3/kg]
            p_gox = PVDW(obj.T_oxtank, obj.m_gox/obj.V_ox_ullage, R_n2o, a_n2o, b_n2o);
        end
        
        function p_oxmanifold = p_oxmanifold(obj)
            % Calculate manifold pressure
            
            % Gravitational acceleration
            g_0 = 9.80665;
            
            % Calculate height difference between liquid level and manifold
            dh_lox = obj.inputs.ox.h_offset_tank + ...
                obj.V_lox/(pi*(obj.inputs.ox.tank_id/2)^2);
            % Calculate acceleration of oxidizer
            if obj.mode.flight_on
                accel_rel = accel_net + g_0;
            else
                accel_rel = g_0;
            end
            % Calculate manifold pressure
            p_oxmanifold = obj.p_oxtank + ...
                accel_rel*obj.N2O_properties.rho_l*dh_lox;
        end
        
        function p_oxtank_press = p_oxtank_press(obj)
            % Calculate oxidizer tank partial pressure of pressurant
            p_oxtank_press = obj.m_oxtank_press*...
                obj.inputs.ox_pressurant.gas_properties.R_specific*...
                obj.T_oxtank/obj.V_ox_ullage;
        end
        
        function p_oxpresstank = p_oxpresstank(obj)
            % Calculate pressure in oxidizer pressurant tank
            m_press_initial = obj.inputs.ox_pressurant.storage_initial_pressure*...
                obj.inputs.ox_pressurant.tank_volume/(...
                obj.inputs.ox_pressurant.gas_properties.R_specific.*...
                obj.inputs.T_amb);
            p_oxpresstank = obj.inputs.ox_pressurant.storage_initial_pressure*...
                (obj.m_oxpresstank/m_press_initial)^...
                obj.inputs.ox_pressurant.gas_properties.gamma;
        end
        
        function T_oxpresstank = T_oxpresstank(obj)
            % Calculate temperature in oxidizer pressurant tank
            T_oxpresstank = obj.inputs.T_amb*...
                (obj.p_oxpresstank/obj.inputs.ox_pressurant.storage_initial_pressure)^...
                ((obj.inputs.ox_pressurant.gas_properties.gamma-1)/...
                obj.inputs.ox_pressurant.gas_properties.gamma);
        end
        
        function V_ox_ullage = V_ox_ullage(obj)
            % Calculate ullage volume in oxidizer tank
            V_ox_ullage = obj.inputs.ox.V_tank - obj.V_lox;
        end
        
        function gamma_ox_ullage = gamma_ox_ullage(obj)
            % Calculate ullage gas ratio of specific heats
            gamma_ox_ullage = (obj.m_gox*obj.N2O_properties.cp_g + ...
                obj.m_oxtank_press*obj.inputs.ox_pressurant.gas_properties.c_p)/...
                (obj.m_gox*obj.N2O_properties.cv_g + ...
                obj.m_oxtank_press*obj.inputs.ox_pressurant.gas_properties.c_v);
        end
        
        function p_cc = p_cc(obj)
            % Calculate pressure in the combustion chamber
            
            R_u = 8.3144621; % Universal gas constant [J/mol*K]
            p_cc = (obj.m_cc/obj.V_cc)*(R_u/obj.M_cc)*obj.T_cc;
        end
    end
    
end

function p = PVDW(T, rho, R, a, b)
%PVDW VanDerWaal's equation of state for pressure
    p = R*T./(1./rho-b)-a*rho.^2;
end