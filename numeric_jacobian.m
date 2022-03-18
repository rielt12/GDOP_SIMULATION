
function jac = numeric_jacobian(f, x, epsilon)
% Calculate Jacobian of function f at given x
% Standard finite difference method
%
% Inputs:
%   f can be a vector of function, but make sure it is a row vector
%   x is where the jacobian is being evaluated, it a row or column vector 
%   epsilon is a very small number

if nargin < 3
    epsilon = 1e-5; 
end

epsilon_inv = 1/epsilon;

nx = length(x); % Dimension of the input x;

f0 = feval(f, x); % caclulate f0, when no perturbation happens

jac = zeros(length(f0), nx);

% Do perturbation
for i = 1 : nx
    xplus = x;
    xplus(i) =  x(i) + epsilon;
    jac(:, i) = (feval(f, xplus) - f0) .* epsilon_inv;
end

    