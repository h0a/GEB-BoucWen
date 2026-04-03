function [x, deltax] = solveLinSys(beam,mesh,loadFactor,x)
% Solving statics problem using a 3-fields mixed formulation and Lagrange multiplier to enforce orthonormal nodal director constraint
% and the Newton-Raphson scheme.
% Input: beam, mesh: struct-variables including information for assembly and solving
%        loadFactor: load factor

% Output: x         = [q; e; s; lagrangeMultipliers] includes the nodal dofs in q, element strains e, element stress s, and Lagrange multipliers
%         deltax    = incremental change of x


% assemble global residual of the mixed formulation
rhs = assembleGlobalResidualStatics(beam, mesh, loadFactor, x);

% assemble global linearized system matrix of the mixed formulation
A = assembleLinSysMatrixStatics(beam, mesh, x, loadFactor, beam.computeAnumerically);

% enforcing essential boundary conditions
[A,rhs] = enforceDirichletBCs(A,rhs,mesh);

% solve
deltax = A \ rhs;

% assemble non-active dofs
deltax = reassembleDirichletBCconstrainedDofs(deltax,mesh);

% update the solution
x = x + deltax;