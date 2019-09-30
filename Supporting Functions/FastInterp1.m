function [ V_q ] = FastInterp1(X, V, X_q)
%FastInterp1 Wrapper for fast use of interp1, working for equally-spaced
%grids.

% Size of grids
n_x = length(X);
dx = (X(end)-X(1))/(n_x-1);

% Get element 1D coordinates
i_x = 1 + (X_q-X(1))/dx;

% Round the coordinates to one high and one low, avoiding going outside the
% available range
[i_x_lo, i_x_hi] = GenerateSamples(i_x, n_x);

% Get high and low values
V_lo = zeros(size(X_q));
V_hi = zeros(size(X_q));
for ii = 1:numel(X_q)
    V_lo(ii) = V(i_x_lo(ii));
    V_hi(ii) = V(i_x_hi(ii));
end

% Interpolate
V_q = (i_x_hi - i_x).*V_lo + (i_x - i_x_lo).*V_hi;
end


function [i_x_lo, i_x_hi] = GenerateSamples(i_x, n_x)
i_x_lo = floor(i_x);
i_x_hi = ceil(i_x);

% Check for integers
i_x_hi(i_x_hi == i_x_lo) = i_x_hi(i_x_hi == i_x_lo) + 1;

% Low range check
i_x_hi(i_x_lo < 1) = 2;
i_x_lo(i_x_lo < 1) = 1;

% High range check
i_x_lo(i_x_hi > n_x) = n_x-1;
i_x_hi(i_x_hi > n_x) = n_x;
end