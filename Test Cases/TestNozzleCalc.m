%% Test NozzleCalc

clear
close all

cd('..\')
addpath(fullfile(pwd, 'Supporting Functions'))

p_up = 100e5; % 100 bar
T_up = 295;
gamma = 1.4;
M = 0.029;
R_u = 8.3144621; % Universal gas constant [J/mol*K]
R_spec = R_u/M;
rho_up = p_up/(R_spec*T_up);
p_back = 0:(p_up/1000):p_up;

d = sqrt(4/pi);

for ii = 1:length(p_back)
    [ ~, ~, ~, ~, m_dot(ii), ~ ] = NozzleCalc(d , d, T_up, p_up, gamma, M, p_back(ii));
    m_dot_choked(ii) = (pi/4*d^2)*p_up/sqrt(T_up)*(gamma/R_spec)^0.5*...
            ((gamma+1)/2)^(-(gamma+1)/(2*(gamma-1)));
    m_dot_inc(ii) = (pi/4*d^2)*sqrt(2*rho_up*(p_up-p_back(ii)));
end

figure
title('NozzleCalc.m Test')
plot(p_back./p_up,m_dot/(pi/4*d^2),p_back./p_up,m_dot_choked/(pi/4*d^2),p_back./p_up,m_dot_inc/(pi/4*d^2))
legend({'NozzleCalc.m','Choked Flow','Incompressible Flow'})
xlabel('P_{down}/P_{up}')
ylabel('Mass Flux [kg/m^2*s]')


cd('.\Test Cases')