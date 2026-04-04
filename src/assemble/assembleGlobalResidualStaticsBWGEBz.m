function rhs = assembleGlobalResidualStaticsBWGEBz(beam, mesh, alpha, x, z)

% current solutions
q = x(1:mesh.ndofs);

% constant external force vector
fext = beam.fext .* alpha;

% residual
resi = zeros(mesh.ndofsTotal,1);

for e = 1:mesh.nelms
    dofsInd = mesh.LM(:,e);
    resi(dofsInd) = resi(dofsInd) + elementBalanceResidual(beam, mesh, e, x);

    dofseInd = mesh.LMstrains(:,e);
    resi(dofseInd) = resi(dofseInd) + elementBWConstitutiveResidual(beam, mesh, e, x, z);

    dofssInd = mesh.LMstresses(:,e);
    resi(dofssInd) = resi(dofssInd) + elementCompatibilityResidual(beam, mesh, e, x);

    dofsLamuInd = mesh.LMlamu(:,e);
    resi(dofsLamuInd) = resi(dofsLamuInd) + elementConstraintResidual(mesh, e, x);
end

rhs = -resi;
rhs(1:mesh.ndofs) = fext + rhs(1:mesh.ndofs);