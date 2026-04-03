function [beam, mesh] = solveStaticFollowingFendGEB(beam,mesh)
    % Solving statics problem of a cantilever subjected to a following end force using a 3-fields mixed formulation 
    % and Lagrange multiplier to enforce orthonormal nodal director constraint and the Newton-Raphson scheme.
    
    % Input: beam, mesh: struct-variables including information for assembly and solving
    % Output:   updated struct-variables: beam, mesh, incl. solutions

    %% ALLOCATION

    % allocation arrays to store solution of primal and dual fields
    beam.q = zeros(mesh.ndofs,beam.numLoadSteps+1);
    beam.q(:,1) = assembleGlobalQ(mesh, beam.nodeVec, beam.nodalD1, beam.nodalD2, beam.nodalD3);
    
    beam.e = zeros(mesh.ndofse,beam.numLoadSteps+1);
    beam.s = zeros(mesh.ndofss,beam.numLoadSteps+1);
    beam.chi = zeros(mesh.numLagrangeMul,beam.numLoadSteps+1);
    
    mesh.num_iters = zeros(beam.numLoadSteps,1);


    %% NEWTON RAPHSON SCHEME

    fprintf('Starting NR iteration...\n');

    % initial guess for the first load step
    x = [beam.q(:,1); beam.e(:,1); beam.s(:,1); beam.chi(:,1)];

    % loop over load steps
    for i = 1:beam.numLoadSteps
        alpha = beam.loadFactors(i);                % load factor    
        NRiter = 0;
    
         while NRiter < beam.maxNumIter
            NRiter = NRiter + 1;
    
            [x, deltax] = solveLinSysStaticFollowingFendGEB(beam,mesh,alpha,x);
    
            if norm(deltax) <= beam.NRtol
                break
            end
        end
    
        % storing NR-converged solutions at the current load step
        mesh.num_iters(i)   = NRiter;
        beam.q(:,i+1)       = x(1:mesh.ndofs);
        beam.e(:,i+1)       = x(mesh.ndofs+1:mesh.ndofs+mesh.ndofse);
        beam.s(:,i+1)       = x(mesh.ndofs+mesh.ndofse+1:mesh.ndofs+mesh.ndofse+mesh.ndofss);
        beam.chi(:,i+1)     = x(mesh.ndofs+mesh.ndofse+mesh.ndofss+1:end);
        
        % stopping computation in case of non-converging NR solution at the current load step
        if NRiter == beam.maxNumIter
            warning('Newton-Raphson scheme does not converge at load step %.1d.\nComputation stopped.\n', i);
            break;
        end
    end
    
    % rescale stress and chi due to scaling factor for conditioning
    beam.s      = beam.s   ./ beam.condScFac;
    beam.chi    = beam.chi ./ beam.condScFac;
    
    % print out number of NR iterations
    [ii, id] = max(mesh.num_iters);
    fprintf('Newton-Raphson scheme finished.\nMax number of NR iterations = %.1d at load step %.1d.\n', ii, id);
end


function [x, deltax] = solveLinSysStaticFollowingFendGEB(beam,mesh,loadFactor,x)
    % Solving the linearized system of equations of a cantilever subjected to a following end force using a 3-fields mixed formulation 
    % and Lagrange multiplier to enforce orthonormal nodal director constraint and the Newton-Raphson scheme.
    
    % Input: beam, mesh: struct-variables including information for assembly and solving
    %        loadFactor: load factor
    %        x:          solution of the previous load step
    
    % Output: x      = [q; e; s; lagrangeMultipliers] includes the nodal dofs in q, element strains e, element stress s, and Lagrange multipliers
    %         deltax = incremental change of x
    
    
    % assemble global residual of the mixed formulation
    rhs = assembleGlobalResidualStatics(beam, mesh, loadFactor, x);
    fextF = assembleFextFollowingF(beam, mesh, loadFactor, x);
    rhs(1:mesh.ndofs) = fextF + rhs(1:mesh.ndofs);
    
    % assemble global linearized system matrix of the mixed formulation
    A = assembleLinSysMatrixStatics(beam, mesh, x, loadFactor, beam.computeAnumerically);
    KfollowingF = aseembleLinTangentKfollowingF(beam, mesh, x, loadFactor, beam.computeAnumerically);
    A(1:mesh.ndofs,1:mesh.ndofs) = -KfollowingF + A(1:mesh.ndofs,1:mesh.ndofs);
    
    % enforcing essential boundary conditions
    [A,rhs] = enforceDirichletBCs(A,rhs,mesh);
    
    % solve
    deltax = A \ rhs;
    
    % assemble non-active dofs
    deltax = reassembleDirichletBCconstrainedDofs(deltax,mesh);
    
    % update the solution
    x = x + deltax;
end