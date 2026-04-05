function [x,NRiter] = solveGEBwithFixedZ(beam, mesh, loadStep, x, z)

alpha = beam.loadFactors(loadStep);                % load factor

NRiter = 0;                                 % iteration count

while NRiter < beam.maxNumIter              % LOOP NR AT CURRENT LOAD STEP
    NRiter = NRiter + 1;

    [x, deltax] = solveLinBWSys(beam,mesh,alpha,x,z);

    if norm(deltax) <= beam.NRtol
        break
    end
end                                         % END LOOP NR AT CURRENT LOAD STEP

% flag to stop computation in case of non-converging NR solution at the current load step
if NRiter == beam.maxNumIter
    warning('Newton-Raphson scheme for GEB does not converge at load step %.1d.\nComputation stopped.\n', loadStep);
    beam.NRflag = 1;
end

end


function [x, deltax] = solveLinBWSys(beam,mesh,loadFactor,x,z)

    % assemble global residual of the mixed formulation
    rhs = assembleGlobalResidualStaticsBWGEBz(beam, mesh, loadFactor, x, z);
    
    % assemble global linearized system matrix of the mixed formulation
    A = assembleLinSysMatrixStaticsBWGEBz(beam, mesh, x, z, loadFactor, beam.computeAnumerically);

    % in case of follower end load or moment
    if beam.nonlinFend
        fextF = assembleFextFollowingF(beam, mesh, loadFactor, x);
        rhs(1:mesh.ndofs) = fextF + rhs(1:mesh.ndofs);
    
        KfollowingF = aseembleLinTangentKfollowingF(beam, mesh, x, loadFactor, beam.computeAnumerically);
        A(1:mesh.ndofs,1:mesh.ndofs) = -KfollowingF + A(1:mesh.ndofs,1:mesh.ndofs);
    end
    
    if beam.nonlinMend
        fextM = assembleFextRollingM(beam, mesh, loadFactor, x);
        rhs(1:mesh.ndofs) = fextM + rhs(1:mesh.ndofs);
    
        KrollingM = aseembleLinTangentKrollingM(beam, mesh, x, loadFactor, beam.computeAnumerically);
        A(1:mesh.ndofs,1:mesh.ndofs) = -KrollingM + A(1:mesh.ndofs,1:mesh.ndofs);
    end
    
    % enforcing essential boundary conditions
    [A,rhs] = enforceDirichletBCs(A,rhs,mesh);
    
    % solve
    deltax = A \ rhs;
    
    % assemble non-active dofs
    deltax = reassembleDirichletBCconstrainedDofs(deltax,mesh);
    
    % update the solution
    x = x + deltax;
end