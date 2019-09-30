function [T] = N2O_FindT(p_vap)
%%N2O_FindT Estimate temperature that gives the specified vapor pressure.

T_min = -90 + 273.15;
T_max = 30 + 273.15;

persistent Pvap

n = 1000;
T_i = linspace(T_min, T_max, n);

if isempty(Pvap)
    Pvap = zeros(size(T_i));

    for ii = 1:n
        prop = N2O_Properties(T_i(ii));
        Pvap(ii) = prop.Pvap;
    end
end

T = interp1(Pvap, T_i, p_vap, 'linear', 'extrap');

end
