function A = assembleLinSysMatrixStaticsBWGEBz(beam, mesh, x, z, alpha, computeAnumerically)


if computeAnumerically
    epsilon = 1e-8;
    A = zeros(mesh.ndofsTotal);
    for i = 1:mesh.ndofsTotal
        x1 = x; x2 = x;
        x1(i) = x1(i) - epsilon;
        x2(i) = x2(i) + epsilon;

        rhs1 = assembleGlobalResidualStaticsBWGEBz(beam, mesh, alpha, x1, z);
        rhs2 = assembleGlobalResidualStaticsBWGEBz(beam, mesh, alpha, x2, z);

        A(:,i) = -(rhs2 - rhs1) ./ (2*epsilon);
    end
    A = sparse(A);
else
    A11 = zeros(mesh.ndofs);    
    A31 = zeros(mesh.ndofss,mesh.ndofs);
    A41 = zeros(mesh.numLagrangeMul,mesh.ndofs);
    A22 = zeros(mesh.ndofse);

    for e = 1:mesh.nelms
        dofsInd = mesh.LM(:,e);
        dofseInd = mesh.LMstrains(:,e) - mesh.ndofs;
        dofssInd = mesh.LMstresses(:,e) - mesh.ndofs - mesh.ndofse;
        dofsLamuInd = mesh.LMlamu(:,e) - mesh.ndofs - mesh.ndofse - mesh.ndofss;

        A11(dofsInd,dofsInd)     = A11(dofsInd,dofsInd)     + elementLinBalanceSysMatrix(beam, mesh, e, x);
        A31(dofssInd,dofsInd)    = A31(dofssInd,dofsInd)    + elementBmatrix(beam, mesh, e, x);
        A41(dofsLamuInd,dofsInd) = A41(dofsLamuInd,dofsInd) + elementHmatrix(mesh, e, x);
        A22(dofseInd,dofseInd)   = A22(dofseInd,dofseInd)   + elementLinBWConstitutiveSysMatrix(beam, mesh, e);
    end
    % employing symmetry
    h = elementArcLength(beam, mesh, 1);    % assuming the same arc length for all elements
    A = sparse( [A11                             sparse(mesh.ndofs,mesh.ndofse)   A31'                      A41'; ...
                 sparse(mesh.ndofse,mesh.ndofs)  A22                             -h .* eye(mesh.ndofse)     sparse(mesh.ndofse,mesh.numLagrangeMul); ...
                 A31                            -h .* eye(mesh.ndofss)            sparse(mesh.ndofss,mesh.ndofss+mesh.numLagrangeMul); ...
                 A41                             sparse(mesh.numLagrangeMul,mesh.ndofse+mesh.ndofss+mesh.numLagrangeMul)] );
end