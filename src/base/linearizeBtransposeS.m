function A = linearizeBtransposeS(s, shapeFuncs, shapeFuncs1stDeriv)
% compute the linearization of (B^T(q) * s) that is a 24x24 matrix A (24 dofs per element)
% input:  s:                    vector s
%         shapeFuncs:           array of the 2 evaluated shape functions L1, L2 (functions associated to qh)
%         shapeFuncs1stDeriv:   array of the evaluated 1st derivative of 2 shape functions L1, L2 (functions associated to qh)
% output: A:                    matrix A of the size 24x24 (24 dofs per element x 24 dofs per element) evaluated at 
%                               the same point as shapeFuncs and shapeFuncs1stDeriv
   
L1 = shapeFuncs(1);
L2 = shapeFuncs(2);
dL1 = shapeFuncs1stDeriv(1);
dL2 = shapeFuncs1stDeriv(2);

D_gamma1 = sparse( blkdiag(dL1 .* eye(3), L1 .* eye(9)) );
D_gamma2 = sparse( blkdiag(dL2 .* eye(3), L2 .* eye(9)) );
D_gamma = [D_gamma1  D_gamma2];

U_gamma = sparse(12,12);
for i = 1:3
    id = 3*i+1:3*(i+1);
    U_gamma(1:3,id) = s(i) .* eye(3);
    U_gamma(id,1:3) = s(i) .* eye(3);
end

D_omega = sparse([ L1 .* eye(12)   L2 .* eye(12); ...
                  dL1 .* eye(12)  dL2 .* eye(12) ]);

U_omega = sparse(24,24);
U_omega(4:6,19:end) = [-s(6).*eye(3) s(5).*eye(3)];
U_omega(19:end,4:6) = [-s(6).*eye(3); s(5).*eye(3)];

U_omega(7:9,16:18) = s(6).*eye(3);
U_omega(7:9,22:24) = -s(4).*eye(3);

U_omega(16:18,7:9) = s(6).*eye(3);
U_omega(22:24,7:9) = -s(4).*eye(3);

U_omega(10:12,16:21) = [-s(5).*eye(3) s(4).*eye(3)];
U_omega(16:21,10:12) = [-s(5).*eye(3); s(4).*eye(3)];

A = D_gamma' * U_gamma * D_gamma + D_omega' * (0.5.*U_omega) * D_omega;