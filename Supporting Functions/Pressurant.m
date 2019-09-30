classdef Pressurant
    %PRESSURANT Represents a propellant pressurant mechanism
    %   Represents an instance of a propellant under pressurization by a
    %   gas.
    
    properties
        % Whether pressurant is to be simulated
        active
        
        % Regulator Pressure Setting (in PSI)
        set_pressure
        % Pressure of Super Charging Tank (in PSI)
        storage_initial_pressure
        % Volume of Super Charging Tank (L)
        tank_volume
        % Regulator CdA (in mm^2)
        flow_CdA
        % Properties of gas
        gas_properties
        
        % Propellant pressurizing ('fuel' or 'oxidizer')
        propellant_type
    end

    methods
        function obj = Pressurant(propellant_type)
            obj.propellant_type = propellant_type;
        end
    end
end
