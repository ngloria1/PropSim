function [crit_flow, mass_flux] = TwoPhaseN2OFlow(T_1, plots_on)
% Calculates the critical mass flow rate using the homogeneous equilibrium 
% model. The mass flow rate may not exceed this, regardless of the back pressure.

% Inputs:
%   T_1: upstream liquid temperature, K
%   plots_on: turn plots on (1) or off (0)

% Outputs:
%   crit_flow:
%       crit_flow.p_up_norm: vector of (upstream pressure / upstream vapor 
%           pressure) for G_crit calculations
%       crit_flow.G_crit: critical oxidizer mass flux at each p_1, kg/m^2*s
%       crit_flow.p_down_norm_crit: critical normalized back pressure 
%           (downstream pressure / upstream pressure) at each p_1, Pa
%   mass_flux:
%       mass_flux.p_up_norm: array of (upstream pressure / upstream vapor 
%           pressure)
%       mass_flux.p_down_norm: array of (back pressure / upstream pressure)
%       mass_flux.G: array of mass flux at corresponding P_1_norm, P_2_norm

if nargin < 2
    plots_on = 0;
end

% Numerical options
n = 200; % density of calculation points along each dimension

% States:
%   1: upstream (stagnation)
%   2: downstream

% saturation properties at upstream temperature
n2o_prop_1_sat = N2O_Properties(T_1);

% Create array of p_1, p_2
p_1_min = 1.0;
p_1_max = 3.0;
p_1 = linspace(p_1_min, p_1_max, n);
p_2 = linspace(0, 1, n);

[P_1_norm,P_2_norm] = meshgrid(p_1,p_2);

% Convert to pressures (from pressure ratio)
P_1 = P_1_norm .* n2o_prop_1_sat.Pvap;
P_2 = P_2_norm .* P_1;

% Upstream liquid enthalpy
% Calculate enthalpy as extension from saturated state
enthalpy_1 = n2o_prop_1_sat.h_l + (P_1 - n2o_prop_1_sat.Pvap)/n2o_prop_1_sat.rho_l;

% Calculate downstream temperature
T_2 = N2O_FindT(P_2);
n2o_prop_2 = N2O_Properties(T_2);

% Calculate gas density
rho_2_g = n2o_prop_2.rho_g;
rho_2_l = n2o_prop_1_sat.rho_l*ones(size(P_2));

enthalpy_2_g = n2o_prop_2.h_g;
enthalpy_2_l = n2o_prop_2.h_l;

entropy_1_l = n2o_prop_1_sat.s_l*ones(size(P_2));
entropy_2_l = n2o_prop_2.s_l;
entropy_2_g = n2o_prop_2.s_g;

% entropy difference: liquid downstream - liquid upstream
entropy_ld_lu = entropy_2_l - entropy_1_l;
% entropy difference: liquid downstream - gas downstream
entropy_ld_gd = entropy_2_l - entropy_2_g;

% calculate mass fraction of gas to conserve entropy
% massfrac: mass fraction of vapor
massfrac = entropy_ld_lu./entropy_ld_gd;
rho_2_l(massfrac < 0) = n2o_prop_1_sat.rho_l;
enthalpy_2_l(massfrac < 0) = n2o_prop_1_sat.h_l...
    + (P_2(massfrac < 0)-n2o_prop_1_sat.Pvap)./n2o_prop_1_sat.rho_l;
massfrac(massfrac < 0) = 0;

% downstream inverse density
rho_2_inv = massfrac.*(1./rho_2_g) + ...
    (1-massfrac).*(1./rho_2_l);

% downstream equivalent density
rho_2_equiv = 1./rho_2_inv;

enthalpy_2 = massfrac.*enthalpy_2_g + ...
    (1-massfrac).*enthalpy_2_l;

% Homogeneous Equilibrium Model
G = rho_2_equiv.*sqrt(2.*(enthalpy_1 - enthalpy_2));

% G_crit for each 
[G_crit, i_crit] = max(G, [], 1);
P_2_crit = zeros(1,size(P_2,2));
for ii = 1:size(P_2,2)
    P_2_crit(ii) = P_2(i_crit(ii),ii);
end

% Create downstream pressure vs. oxidizer mass flux profile
G_out = repmat(G_crit, [size(P_2,1), 1]);
% P_2_crit expanded to size of P_2
P_2_crit_exp = repmat(P_2_crit, [size(P_2,1), 1]);
G_out(P_2 > P_2_crit_exp) = G(P_2 > P_2_crit_exp);

% Calculate incompressible mass flux
G_inc = sqrt(2*(P_1-P_2)*n2o_prop_1_sat.rho_l);


% Plotting
n_plots = 3;
if plots_on
    for plot_index = round(linspace(1,size(P_2,2),n_plots))
        P_1_norm_index = P_1_norm(:,plot_index);
        P_1_norm_index = P_1_norm_index(~isnan(P_1_norm_index));
        P_1_norm_index = mean(P_1_norm_index);

        figure
        subplot(2,2,1)
        plot(P_2_norm(:,plot_index), G(:,plot_index),...
            P_2_norm(:,plot_index), G_out(:,plot_index),  ...
            P_2_norm(:,plot_index), G_inc(:,plot_index))
        xlabel('Normalized Back Pressure (p_{downstream}/p_{upstream})')
        ylabel('Ox Mass Flux [kg/m^2*s]')
        xlim([min(P_2_norm(:,plot_index)),max(P_2_norm(:,plot_index))])
        legend({'HEM - Calc.','HEM - Physical','Incompressible'})
        subplot(2,2,2)
        plot(P_2_norm(:,plot_index), G_out(:,plot_index))
        xlabel('Normalized Back Pressure (p_{downstream}/p_{upstream})')
        ylabel('Ox Mass Flux [kg/m^2*s]')
        xlim([min(P_2_norm(:,plot_index)),max(P_2_norm(:,plot_index))])
        subplot(2,2,3)
        plot(P_2_norm(:,plot_index), massfrac(:,plot_index))
        xlabel('Normalized Back Pressure (p_{downstream}/p_{upstream})')
        ylabel('Vapor Mass Fraction')
        subplot(2,2,4)
        plot(P_2_norm(:,plot_index), rho_2_equiv(:,plot_index), ...
            P_2_norm(:,plot_index), rho_2_g(:,plot_index), ...
            P_2_norm(:,plot_index), rho_2_l(:,plot_index))
        xlabel('Normalized Back Pressure (p_{downstream}/p_{upstream})')
        ylabel('Density [kg/m^3]')
        legend({'Effective','Gas','Liquid'})
        suptitle(sprintf('p_{upstream}/p_{vapor} = %.3g', P_1_norm_index))
    end
end

% Package for output
crit_flow.p_up_norm = p_1;
crit_flow.G_crit = G_crit;
crit_flow.p_down_norm_crit = P_2_crit./(p_1*n2o_prop_1_sat.Pvap);
mass_flux.p_up_norm = P_1_norm;
mass_flux.p_down_norm = P_2_norm;
mass_flux.G = G_out;
end
