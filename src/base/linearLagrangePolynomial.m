function L = linearLagrangePolynomial(a,b)

% computing linear Lagrange polynomial and their derivatives on an arbitrary interval [a,b]
% based on the definition on [0,1] - wikipedia

% Input: interval [a,b]
% Output: L{i} is array of the (i-1)-th order of functions {L1 L2}

L = cell(2,1);

% p(t) = L1 * p1 + L2 * p2

% o-th derivative
L1 = @(x)  (x-b) / (a-b);
L2 = @(x)  (x-a) / (b-a);

L{1} = {L1, L2};


% 1st derivative
dL1 = @(x) 1 / (a-b); 
dL2 = @(x) 1 / (b-a);

L{2} = {dL1, dL2};
