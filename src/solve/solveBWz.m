function [z2,BWiter] = solveBWz(beam, mesh, loadStep, x1, x2, z1, z2)
% solving for z of each strain dof using the Newton-Raphson method when
% using the Bouc-Wen hysteretic model

% Input:  beam:      struct-variable with all input parameters
%         loadStep:  id od the current load step
%         x1:        the (i)-th   structural solution
%         x2:        the (i+1)-th structural solution
%         z1:        the (i)-th   z-solution for each strain dof
%         z2:        the (i+1)-th z-solution for each strain dof
% Output: z2:        converged (i+1)-th z-solution for each strain dof
%         maxBWiter: max number of NR iteration over all strain dofs needed for this solving procedure


% allocation BWiter to get max BWiter for output
BWiter = [];

% solving for each z-component corresponding to each strain dof, given that
% alpha4Sigma_i is not 1 (if alpha4Sigma_i = 1, this i-th z-component is
% not relevant)
for e = 1:mesh.nelms                          % LOOP OVER ELEMENTS
    % structural strain field of the current element
    e1 = x1(mesh.LMstrains(:,e));
    e2 = x2(mesh.LMstrains(:,e));

    % corresponding z of the current element
    z1e = z1(mesh.LMstrains(:,e)-mesh.ndofs);
    z2e = z2(mesh.LMstrains(:,e)-mesh.ndofs);

    % solving for relevant z-component
    for i = 1:numel(z1e)                      % LOOP OVER Z-COMPONENTS OF THE CURRENT ELEMENT
        if (beam.alpha4Sigma(i) == 1) == false
            [zi,iter] = NewtonRaphsonScheme4z(beam, e1(i), e2(i), z1e(i), z2e(i));

            % overwrite with converged z and store number of NR iterations
            z2e(i) = zi;
            BWiter = [BWiter;iter];

            % flag to stop computation in case of non-converging NR solution
            if iter == beam.maxNumBWIter
                warning('Newton-Raphson scheme for zBW does not converge at load step %.1d.\nComputation stopped.\n', loadStep);
                beam.BWflag = 1; break
            end
        end
    end                                       % END LOOP OVER Z-COMPONENTS OF THE CURRENT ELEMENT

    % overwrite with converged z
    z2(mesh.LMstrains(:,e)-mesh.ndofs) = z2e;

    if beam.BWflag == 1
        break
    end
end                                           % END LOOP OVER ELEMENTS

% get max number of iterations over all strain dofs
if isempty(BWiter)
    BWiter = 0; % Initialize BWiter to zero if no iterations were recorded
else
    BWiter = max(BWiter);
end
end




function [z2,BWiter] = NewtonRaphsonScheme4z(beam, e1, e2, z1, z2)

    BWiter = 0;                                   % iteration count
    
    while BWiter < beam.maxNumBWIter              % LOOP NR FOR CURRENT Z-COMPONENT
        BWiter = BWiter + 1;

        Jz = Jacobian_2 (beam.BWdata, e2, e1, z1, z2);
        hz = Residuum_2 (beam.BWdata, e2, e1, z1, z2);
        deltaz = -hz / Jz;

        % update z
        z2 = z2 + deltaz;

        % check convergence
        if abs(deltaz) <= beam.BWNRtol
            break
        end
    end                                           % END LOOP NR FOR CURRENT Z-COMPONENT
end
