function [beam, mesh] = solveStaticGEB(beam, mesh)

% Solving statics problem of GEB, using a 3-fields mixed formulation and Lagrange multiplier method to enforce orthonormal nodal director constraint
% and without using the nullspace approach.
% Solver algorithm is based on the Newton-Raphson scheme for the nonlinear formulation of GEB.

% Input:    beam, mesh: struct-variables including information for assembly and solving
% Output:   updated struct-variables: beam, mesh, incl. solutions



%% ALLOCATION

% allocation arrays to store solution of primal and dual fields
beam.q = zeros(mesh.ndofs,beam.numLoadSteps+1);
beam.q(:,1) = assembleGlobalQ(mesh, beam.nodeVec, beam.nodalD1, beam.nodalD2, beam.nodalD3);

beam.e = zeros(mesh.ndofse,beam.numLoadSteps+1);
beam.s = zeros(mesh.ndofss,beam.numLoadSteps+1);
beam.chi = zeros(mesh.numLagrangeMul,beam.numLoadSteps+1);

mesh.num_iters = zeros(beam.numLoadSteps,1);

% external force vector (constant direction)
beam.fext = assembleFextNodal(beam, mesh) + assembleFextBody(beam, mesh);



%% NEWTON RAPHSON SCHEME

fprintf('Starting NR iteration...\n');

% initial guess for the 1st load step
x = [beam.q(:,1); beam.e(:,1); beam.s(:,1); beam.chi(:,1)];

% loop over load steps
for i = 1:beam.numLoadSteps                     % LOOP OVER LOAD STEPS
    alpha = beam.loadFactors(i);                % load factor 
    
    % solving with Newton Raphson
    [x,beam,mesh] = NewtonRaphsonScheme(beam,mesh,x,alpha,i);
    
    % stopping computation in case of non-converging NR solution at the current load step
    if beam.NRflag == 1
        break;
    end
end                                             % END LOOP OVER LOAD STEPS

% rescale stress and chi due to scaling factor for conditioning
beam.s      = beam.s   ./ beam.condScFac;
beam.chi    = beam.chi ./ beam.condScFac;

% print out number of NR iterations
[ii, id] = max(mesh.num_iters);
fprintf('Newton-Raphson scheme finished.\nMax number of NR iterations = %.1d at load step %.1d.\n', ii, id);
