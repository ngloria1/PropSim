%% Plot N2O Properties

close all
clear

addpath(fullfile('..', 'Supporting Functions'))

% Plot Nitrous Properties over Temperature Range
n = 1000;
T = linspace(-90,30,n) + 273.15;

cp_l = zeros(n,1);
cv_l = zeros(n,1);
cp_g = zeros(n,1); 
cv_g = zeros(n,1);
Pvap = zeros(n,1);
dH_vap = zeros(n,1);
dE_vap = zeros(n,1);
rho_g = zeros(n,1);
rho_l = zeros(n,1);

for i = 1:n
    Prop = N2O_Properties(T(i));
    cp_l(i) = Prop.cp_l;
    cv_l(i) = Prop.cv_l;
    cp_g(i) = Prop.cp_g;
    cv_g(i) = Prop.cv_g;
    Pvap(i) = Prop.Pvap;
    dH_vap(i) = Prop.deltaH_vap;
    dE_vap(i) = Prop.deltaE_vap;
    rho_g(i) = Prop.rho_g;
    rho_l(i) = Prop.rho_l;
end
gamma = cp_g./cv_g;

% Density
figure
plot(T - 273.15,rho_l, T - 273.15, rho_g)
xlabel('Temperature [C]')
ylabel('Density [kg/m^3]')
title('Nitrous Oxide Density')

% Specific Heats / dH_vap
figure
subplot(2,2,1)
plot(T - 273.15,cv_l,T - 273.15,cp_l)
title('Liquid Specific Heat')
xlabel('Temperature [deg C]')
ylabel('Specific Heat [J/kg*K]')
legend({'Isochoric','Isobaric'})
subplot(2,2,2)
plot(T - 273.15,cv_g,T - 273.15,cp_g)
title('Gas Specific Heats')
xlabel('Temperature [deg C]')
ylabel('Specific Heat [J/kg*K]')
legend({'Isochoric','Isobaric'})
subplot(2,2,3)
plot(T - 273.15,gamma)
title('Gas Ratio of Specific Heats')
xlabel('Temperature [deg C]')
ylabel('\gamma []')
subplot(2,2,4)
plot(T - 273.15,dH_vap,T - 273.15,dE_vap)
title('Enthalpy of Vaporization')
xlabel('Temperature [deg C]')
ylabel('Enthalpy of Vaporization [J/kg*K]')
legend({'dH','dE'})

% Vapor Pressure
figure
plot(T - 273.15,Pvap/1e6)
title('Vapor Pressure')
xlabel('Temperature [deg C]')
ylabel('Vapor Pressure [MPa]')

