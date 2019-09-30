function [ ischoked, p_exit, p_shock_crit, u_exit, m_dot, M_e ] = ...
    NozzleCalc( d_throat, d_exit, T0, p0, gamma, M, p_back)
%NozzleCalc Calculates the exit conditions and mass flow rate for flow 
%through a nozzle. Assumes adiabatic flow. 
%   INPUTS:
%       - d_throat: throat diameter [m]
%       - d_exit: nozzle exit diameter [m]
%       - T0: stagnation temperatre [K]
%       - p0: stagnation pressure at nozzle inlet [Pa]
%       - gamma: ratio of specific heats []
%       - M: mean molecular mass [kg/mol]
%       - p_back: exit back pressure [Pa]
%   OUTPUS:
%       - ischoked: 0 for not choked, 1 for choked
%       - p_exit: static pressure at exit [Pa]
%       - p_shock_crit: the back pressure above which a shock occurs [Pa]
%       - u_exit: axial velocity at exit [m/s]
%       - m_dot: mass flow rate through nozzle [kg/s]
%       - M_e: exit mach number []

A_throat = pi/4*d_throat^2;
A_exit = pi/4*d_exit^2;
E = A_exit/A_throat; % Expansion ratio []

R_u = 8.3144621; % Universal gas constant [J/mol*K]
R_spec = R_u/M;

rho0 = p0/(R_spec*T0);

%% Determine Regime of Nozzle Flow
% Calculate Shock Critical Pressure
p_shock_crit = FindShockPCrit(A_throat, A_exit, gamma, p0);
if p_back > FindChokedPCrit(gamma, p0) 
    % All subsonic
    ischoked = 0;
    p_exit = p_back;
    p_rat = p_back/p0;
    [M_e,T_rat,~,rho_rat,~] = flowisentropic(gamma, p_rat, 'pres');
    T = T_rat*T0;
    a = sqrt(gamma*R_spec*T);
    u_exit = M_e*a;
    rho = rho0*rho_rat;
    m_dot = A_throat*u_exit*rho;
elseif p_back < p_shock_crit 
    % Choked flow, supersonic at exit
    ischoked = 1;
    [M_e,T_e_rat,p_e_rat,rho_e_rat,~] = flowisentropic(gamma, E, 'sup');
    p_exit = p0*p_e_rat;
    T_e = T0*T_e_rat;
    rho_e = rho0*rho_e_rat;
    a_e = sqrt(gamma*R_spec*T_e);
    u_exit = a_e*M_e;
    m_dot = u_exit*rho_e*A_exit;
else
    % Choked flow, normal shock in diverging section
    ischoked = 1;
    options = optimoptions('fsolve', 'TolX', 1e-3, 'Display', 'off');
    E_shock = fsolve( @SolveNormShock, 0.5*E, options);
    [p_exit, u_exit, m_dot, M_e] = NormalShockCalc( A_throat, E_shock*A_throat, ...
    A_exit, gamma, R_spec, p0, T0);
end

function [value] = SolveNormShock( E_shock )
    % Compute pressure difference between exit and ambient for shock in
    % diverging part of nozzle. 
    % Conditions:
    % - 0: period before shock
    % - 1: right before shock
    % - 2: right after shock
    % - 3: at nozzle exit
    % - B: period after shock

    [p_exit_guess, ~, ~, ~] = NormalShockCalc( A_throat, ...
        E_shock*A_throat, A_exit, gamma, R_spec, p0, T0);
    value = p_exit_guess - p_back;
end

end

function [p_exit] = FindChokedPCrit(gamma, p0)
    % Find cut-off ambient pressure for choked flow
    p_rat = (2/(gamma+1))^(gamma/(gamma-1));
    p_exit = p_rat*p0;
end

function [p_exit] = FindShockPCrit(A_throat, A_exit, gamma, p0)
    % Find cut-off ambient pressure for shock in nozzle
    A_star = A_throat;
    A_rat = A_exit/A_star;
    % Limiting condition: normal shock at exit plane
    [M_e,~,p_rat_isen,~,~] = flowisentropic(gamma, A_rat, 'sup');
    p_exit_before = p_rat_isen*p0;
    [~,~,p_rat_shock,~,~,~,~] = flownormalshock(gamma, M_e, 'mach');
    p_exit = p_rat_shock*p_exit_before;
end

function [p_exit, u_exit, m_dot, M_e] = NormalShockCalc( A_throat, A_shock, ...
    A_exit, gamma, R_spec, p0, T0)
% Compute mass flow rate and exit pressure for normal shock at some point
% in diverging section of nozzle based on the area ratio of shock location.
    % E_shock: A_shock/A_throat
    % Conditions:
    % - A: period before shock
    % - 1: right before shock
    % - 2: right after shock
    % - e: at nozzle exit
    % - B: period after shock
    A_shocktoAstar_rat = A_shock/A_throat;
    A_exittoAstar_rat = A_exit/A_throat;
    
    % Calculate conditions at 1
    [M_1,T_1_rat,p_1_rat,~,~] = flowisentropic(gamma, A_shocktoAstar_rat, 'sup');
    T_1 = T0*T_1_rat;
    p_1 = p0*p_1_rat;
    % Calculate conditions at 2
    [~,p_2_rat,T_2_rat,~,M_2,p0B_rat,~] = flownormalshock(gamma, M_1, 'mach');
    p_2 = p_1*p_2_rat;
    T_2 = T_1*T_2_rat;
    p0B = p0B_rat*p0;
    % Calculate E for new A_star
    [~,~,~,~,A_2toBstar_rat] = flowisentropic(gamma, M_2, 'mach');
    A_exittoBstar_rat = A_exittoAstar_rat/A_shocktoAstar_rat*A_2toBstar_rat; % A_exit/A_star,B
    % Calculate conditions at 3
    [M_e,T_e_rat,p_e_rat,~,~] = flowisentropic(gamma, A_exittoBstar_rat, 'sub');
    T_e = T0*T_e_rat;
    p_exit = p0B*p_e_rat;
    rho_e = p_exit/(R_spec*T_e);
    a_e = sqrt(gamma*R_spec*T_e);
    u_exit = M_e*a_e;
    m_dot = u_exit*rho_e*A_throat;
end