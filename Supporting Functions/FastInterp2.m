function [ V_q ] = FastInterp2(X, Y, V, X_q, Y_q)
%FastInterp2 Wrapper for fast use of interp2, working for equally-spaced
%grids.

% Size of grids
n_x = size(X,2);
n_y = size(Y,1);

% Get 1D values for each dimension
X_i = squeeze(X(1,1:n_x,1));
Y_i = squeeze(Y(1:n_y,1,1));

% Get element 1D coordinates
i_x = FastInterp1(X_i,1:n_x,X_q);
i_y = FastInterp1(Y_i,1:n_y,Y_q);

% Round the coordinates to one high and one low, avoiding going outside the
% available range
[i_x_lo, i_x_hi] = GenerateSamples(i_x, n_x);
[i_y_lo, i_y_hi] = GenerateSamples(i_y, n_y);

% Reduce matrices to neigboring units
V_red = V(i_y_lo:i_y_hi,i_x_lo:i_x_hi);
X_red = X(i_y_lo:i_y_hi,i_x_lo:i_x_hi);
Y_red = Y(i_y_lo:i_y_hi,i_x_lo:i_x_hi);

% Find the 
V_q = interp2(X_red, Y_red, V_red, X_q, Y_q);
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