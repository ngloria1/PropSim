function [ properties ] = N2O_Properties( T )
% N2O_Properties( T )
% Calculates properties of nitrous 
% oxide given a temperature in K
% WARNING: if temperature input is outside of -90 to 30 C, properties will
% be generated for boundary (-90 C or 30 C)
% 
% Returns properties (structure)
% properties.Pvap in Pa
% properties.rho_l in kg/m^3
% properties.rho_g in kg/m^3
% properties.deltaH_vap in J/kg
% properties.cp_l in J/(kg*K)
% properties.cv_l in J/(kg*K)
% properties.cp_g in J/(kg*K)
% properties.cv_g in J/(kg*K)
% properties.h_l in J/(kg*K)
% properties.h_g in J/(kg*K)
% properties.s_l in J/(kg*K)
% properties.s_g in J/(kg*K)
% properties.mu_l in N*s/(m^2)
% properties.mu_g in N*s/(m^2)

persistent NIST_data

R_u = 8.3144621; % Universal gas constant [J/mol*K]
M_n2o = 0.044013; % Molecular mass of nitrous oxide [kg/mol]
R_n2o_0 = R_u/M_n2o; % Specific gas constant of nitrous oxide [J/kg*K]

% Range-check temperature
T(T < (-90 + 273.15)) = -90 + 273.150001;
T(T > (30 + 273.15)) = 30 + 273.150001;

Tcrit = 309.57; %K
Pcrit = 7251;   %kPa
rhocrit = 452; %kg/m^3
%possibly add critical compressibility factor "Z"
Tr = T/Tcrit;

% Calculate vapor pressure, valid -90 to 36C
b1 = -6.71893;
b2 = 1.3596;
b3 = -1.3779;
b4 = -4.051;

properties.Pvap = exp((1./Tr).*(b1*(1-Tr) + b2*(1-Tr).^(3/2) + b3*(1-Tr).^(5/2) + b4*(1-Tr).^5))*Pcrit;
properties.Pvap = properties.Pvap*1000;

% Calculate Density of Liquid, valid -90C to 36C
b1 = 1.72328;
b2 = -0.83950;
b3 = 0.51060;
b4 = -0.10412;

properties.rho_l = exp(b1*(1-Tr).^(1/3) + b2*(1-Tr).^(2/3) + b3*(1-Tr) + b4*(1-Tr).^(4/3))*rhocrit;

% Calculate Density of Gas, valid -90C to 36C
b1 = -1.00900;
b2 = -6.28792;
b3 = 7.50332;
b4 = -7.90463;
b5 = 0.629427;
Trinv = 1./Tr;

properties.rho_g = exp(b1*(Trinv-1).^(1/3) + b2*(Trinv-1).^(2/3) + b3*(Trinv-1) + b4*(Trinv-1).^(4/3) + b5*(Trinv-1).^(5/3))*rhocrit;

% Calculate dynamic viscosity of saturated liquid, valid from -90C to 30C
b1 = 1.6089;
b2 = 2.0439;
b3 = 5.24;
b4 = 0.0293423;
theta = (Tcrit-b3)./(T-b3);

properties.mu_l = b4*exp(b1*(theta-1).^(1/3) + b2*(theta-1).^(4/3));

% Calculate dynamic viscosity of saturated vapor, valid from -90C to 30C
b1 = 3.3281;
b2 = -1.18237;
b3 = -0.055155;
Trinv = 1./Tr;

properties.mu_g = exp(b1 + b2*(Trinv-1).^(1/3) + b3*(Trinv-1).^(4/3));

% Find Specific Enthalpy
if isempty(NIST_data)
    data = tdfread('N2O_Properties.cgi.txt');
    NIST_data.T = data.Temperature_0x28K0x29;
    NIST_data.h_liq = data.Enthalpy_0x28l0x2C_kJ0x2Fkg0x29;
    NIST_data.h_gas = data.Enthalpy_0x28v0x2C_kJ0x2Fkg0x29;
    NIST_data.e_liq = data.Internal_Energy_0x28l0x2C_kJ0x2Fkg0x29;
    NIST_data.e_gas = data.Internal_Energy_0x28v0x2C_kJ0x2Fkg0x29;
    NIST_data.cv_l = data.Cv_0x28l0x2C_J0x2Fg0x2AK0x29; % Cv for liquid
    NIST_data.cv_g = data.Cv_0x28v0x2C_J0x2Fg0x2AK0x29; % Cv for gas
    NIST_data.cp_l = data.Cp_0x28l0x2C_J0x2Fg0x2AK0x29; % Cp for liquid
    NIST_data.cp_g = data.Cp_0x28v0x2C_J0x2Fg0x2AK0x29; % Cp for gas
    NIST_data.s_l = data.Entropy_0x28l0x2C_J0x2Fg0x2AK0x29;
    NIST_data.s_g = data.Entropy_0x28v0x2C_J0x2Fg0x2AK0x29;
end
% Gas Specific Enthalpy
properties.h_l = FastInterp1(NIST_data.T, NIST_data.h_liq, T)*1e3; % J/kg
properties.h_g = FastInterp1(NIST_data.T, NIST_data.h_gas, T)*1e3; % J/kg
properties.e_l = FastInterp1(NIST_data.T, NIST_data.e_liq, T)*1e3; % J/kg
properties.e_g = FastInterp1(NIST_data.T, NIST_data.e_gas, T)*1e3; % J/kg
properties.deltaH_vap = properties.h_g-properties.h_l;
properties.deltaE_vap = properties.e_g-properties.e_l;

% Specific Heat at Constant Volume
properties.cv_l = FastInterp1(NIST_data.T, NIST_data.cv_l, T)*1e3;
properties.cv_g = FastInterp1(NIST_data.T, NIST_data.cv_g, T)*1e3;

% Specific Heat at Constant Pressure
properties.cp_l = FastInterp1(NIST_data.T, NIST_data.cp_l, T)*1e3;
properties.cp_g = FastInterp1(NIST_data.T, NIST_data.cp_g, T)*1e3;

%Specific Entropy

properties.s_l = FastInterp1(NIST_data.T, NIST_data.s_l, T)*1e3;
properties.s_g = FastInterp1(NIST_data.T, NIST_data.s_g, T)*1e3;

%% Convert Properties to Standard Units
properties.mu_l = properties.mu_l*10^-3; % mN*s/(m^2) -> N*s/m^2
properties.mu_g = properties.mu_g*10^-6; % uN*s/(m^2) -> N*s/m^2

end
