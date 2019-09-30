function [G, G_crit, p_down_crit] = FindG(T_up, p_up, p_down)
% Find critical mass flux for injector based on stored data from two-phase 
% flow calculator function
% Inputs:
%   T_up: upstream oxidizer temperature [K]
%   p_up: upstream pressure [Pa]
%   p_down: downstream pressure [Pa]
% Outputs:
%   G: mass flux [kg/m^2*s]
%   G_crit: critical mass flux for specified upstream pressure [kg/m^2*s]
%   p_down_crit: critical downstream pressure for specified upstream 
%       pressure [Pa]

persistent crit_flow T_ii p_1_norm_ii p_2_norm_ii G_i ranges

n2o_prop = N2O_Properties(T_up);
p_vap = n2o_prop.Pvap;

if isempty(G_i)
    %% Initialize array of injector critical mass flux values
    
    % Initialize Grid
    ranges.T_range = [-90, 30] + 273.15; 
    ranges.p_1_range = [1, 3]; 
    ranges.p_2_range = [0, 1]; 
    
    % grid size
    n = 200;
    
    T_i = linspace(ranges.T_range(1),ranges.T_range(2),n);
    p_1_norm_i = linspace(ranges.p_1_range(1),ranges.p_1_range(2),n);
    p_2_norm_i = linspace(ranges.p_2_range(1),ranges.p_2_range(2),n); 

    [T_ii, p_1_norm_ii, p_2_norm_ii] = meshgrid(T_i, p_1_norm_i, p_2_norm_i);
    [crit_flow.T_ii, crit_flow.p_1_norm_ii] = meshgrid(T_i, p_1_norm_i);

    fprintf('Calculating Oxidizer Mass Flux Array...\n');

    G_i = zeros(size(T_ii));
    crit_flow.G_crit = zeros(size(crit_flow.T_ii));
    for ii = 1:length(T_i)
        [crit_flow_i, mass_flux] = TwoPhaseN2OFlow(T_i(ii));
        crit_flow.G_crit(:,ii) = interp1(crit_flow_i.p_up_norm,...
            crit_flow_i.G_crit,p_1_norm_i);
        crit_flow.p_down_norm_crit(:,ii) = interp1(crit_flow_i.p_up_norm,...
            crit_flow_i.p_down_norm_crit,p_1_norm_i);
        
        G_i(:,ii,:) = interp2(mass_flux.p_up_norm, ...
            mass_flux.p_down_norm, mass_flux.G, p_1_norm_ii(:,ii,:), ...
            p_2_norm_ii(:,ii,:));
    end

end
% Interpolation error catching
if (T_up > ranges.T_range(2))||(T_up < ranges.T_range(1))
    warning('T_up = %f out of range: [%f, %f]', T_up, ...
        ranges.T_range(1), ranges.T_range(2))
    T_up = max(min(T_up,ranges.T_range(2)),ranges.T_range(1));
elseif (p_up/p_vap > ranges.p_1_range(2))||...
        (p_up/p_vap < ranges.p_1_range(1))
    warning('p_up/p_vap = %f out of range: [%f, %f]', p_up/p_vap, ...
        ranges.p_1_range(1), ranges.p_1_range(2))
    p_up = max(min(p_up,p_vap*ranges.p_1_range(2)),p_vap*ranges.p_1_range(1));
elseif (p_down/p_up > ranges.p_2_range(2))||...
        (p_down/p_up < ranges.p_2_range(1))
    warning('p_down/p_up = %f out of range: [%f, %f]', p_down/p_up, ...
        ranges.p_2_range(1), ranges.p_2_range(2))
    p_down = max(min(p_down,p_up*ranges.p_2_range(2)),p_up*ranges.p_2_range(1));
end

% Interpolation
 G = FastInterp3(T_ii, p_1_norm_ii, p_2_norm_ii, G_i, T_up, p_up/p_vap, p_down/p_up);
 G_crit = FastInterp2(crit_flow.T_ii,crit_flow.p_1_norm_ii,...
     crit_flow.G_crit,T_up,p_up/p_vap);
 p_down_crit = FastInterp2(crit_flow.T_ii,crit_flow.p_1_norm_ii,...
     crit_flow.p_down_norm_crit,T_up,p_up/p_vap)*p_up;
end
