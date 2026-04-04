function [beam, mesh] = solveStaticBWGEB(beam, mesh)

% Solving statics problem of GEB, using a 3-fields mixed formulation and Lagrange multiplier method to enforce orthonormal nodal director constraint
% and the Bouc-Wen hysteretic model, without using the nullspace approach.
% Solver algorithm is based on the Newton-Raphson scheme for the nonlinear formulation of GEB and Uzawa method.

% Input:    beam, mesh: struct-variables including information for assembly and solving
% Output:   updated struct-variables: beam, mesh, incl. solutions



%% ALLOCATION

% allocation arrays to store solution of primal and dual fields
beam.q = zeros(mesh.ndofs,beam.numLoadSteps+1);
beam.q(:,1) = assembleGlobalQ(mesh, beam.nodeVec, beam.nodalD1, beam.nodalD2, beam.nodalD3);

beam.e = zeros(mesh.ndofse,beam.numLoadSteps+1);
beam.e(:,1) = beam.e0;      % initial condition

beam.s = zeros(mesh.ndofss,beam.numLoadSteps+1);
beam.chi = zeros(mesh.numLagrangeMul,beam.numLoadSteps+1);

beam.zBW = zeros(mesh.ndofse,beam.numLoadSteps+1);
beam.zBW(:,1) = beam.z0;    % initial condition

mesh.num_iters      = zeros(beam.numLoadSteps,beam.maxUzawaIter);
mesh.num_iters4z    = zeros(beam.numLoadSteps,beam.maxUzawaIter);
mesh.num_itersUzawa = zeros(beam.numLoadSteps,1);
mesh.errUzawa       = zeros(beam.maxUzawaIter,beam.numLoadSteps,3);

% external force vector (constant direction)
beam.fext = assembleFextNodal(beam, mesh) + assembleFextBody(beam, mesh);



%% SOLVER

fprintf ('-------------------------------------------------------------\n')
fprintf ('Static analysis of GEB with the Bouc-Wen hysteretic model \n')
fprintf ('-------------------------------------------------------------\n')
fprintf ('Time Step h = %8.6f \n', beam.timeStep)
fprintf ('Starting computations over %i time/load steps...\n',beam.numLoadSteps);
fprintf ('-------------------------------------------------------------\n')


for i = 1:beam.numLoadSteps                     % LOOP OVER LOAD STEPS    

    % get the (i)-solution
    x1 = [beam.q(:,i); beam.e(:,i); beam.s(:,i); beam.chi(:,i)];
    z1 = beam.zBW(:,i);

    % initial guess of the (i+1)-solution
    x2 = x1;
    z2 = z1;

    % counting Uzawa iterations
    UzawaIter = 1;

    % Uzawa loop for the (i+1)-solution
    while UzawaIter < beam.maxUzawaIter        % LOOP UZAWA SCHEME
        % get the (i+1)-solution of the previous Uzawa iteration
        xPre = x2;
        zPre = z2;

        % solve GEB
        [x2,NRiter] = solveGEBwithFixedZ(beam, mesh, i, x2, z2);
        mesh.num_iters(i,UzawaIter)   = NRiter;
        if beam.NRflag == 1
            break;
        end

        % solve zBW
        [z2,BWiter] = solveBWz(beam, mesh, i, x1, x2, z1, z2);
        mesh.num_iters4z(i,UzawaIter)   = BWiter;
        if beam.BWflag == 1
            break;
        end

        % relaxation for zBW
        z2 = relaxBWz(beam, mesh, zPre, z2);

        % compute difference in x and zBW
        err = max( norm(x2-xPre), norm(z2-zPre) );
        mesh.errUzawa(UzawaIter,i,:) = [norm(x2-xPre), norm(z2-zPre), err];

        % check convergence of the Uzawa scheme
        if err <= beam.UzawaTol
            break;
        else
            UzawaIter = UzawaIter + 1;
        end
    end                                         % END LOOP UZAWA SCHEME

    % storing converged solutions at the current load step
    mesh.num_itersUzawa(i) = UzawaIter;

    beam.q(:,i+1)       = x2(1:mesh.ndofs);
    beam.e(:,i+1)       = x2(mesh.ndofs+1:mesh.ndofs+mesh.ndofse);
    beam.s(:,i+1)       = x2(mesh.ndofs+mesh.ndofse+1:mesh.ndofs+mesh.ndofse+mesh.ndofss);
    beam.chi(:,i+1)     = x2(mesh.ndofs+mesh.ndofse+mesh.ndofss+1:end);
    
    beam.zBW(:,i+1)     = z2;

    % stop computations if any scheme does not converge
    if (beam.NRflag == 1) || (beam.BWflag == 1) || (UzawaIter == beam.maxUzawaIter)

        if UzawaIter == beam.maxUzawaIter
            warning('The Uzawa scheme does not converge at load step %.1d.\nComputation stopped.\n', i);
        end

        break;      % break load step loop
    end
end                                             % END LOOP OVER LOAD STEPS



%% POSTPROCESSING

% rescale stress and chi due to scaling factor for conditioning
beam.s      = beam.s   ./ beam.condScFac;
beam.chi    = beam.chi ./ beam.condScFac;

% print out max number of iterations
[ii, id] = max(mesh.num_itersUzawa);
fprintf('Computation finished.\nMax number of Uzawa iterations = %.1d at load step %.1d.\n', ii, id);

[maxNRiter, linIdx] = max(mesh.num_iters(:));
[loadStep, UzIter] = ind2sub(size(mesh.num_iters), linIdx);
fprintf('Max number of GEB-NR iterations = %.1d at load step %.1d and %.1d-th Uzawa iteration.\n', maxNRiter, loadStep, UzIter);

[maxBWiter, linIdx] = max(mesh.num_iters4z(:));
[loadStep, UzIter] = ind2sub(size(mesh.num_iters4z), linIdx);
fprintf('Max number of Bouc-Wen-NR iterations = %.1d at load step %.1d and %.1d-th Uzawa iteration.\n', maxBWiter, loadStep, UzIter);
