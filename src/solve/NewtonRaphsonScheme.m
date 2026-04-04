function [x,beam,mesh] = NewtonRaphsonScheme(beam,mesh,x,loadFac,loadStep)

alpha = loadFac;
i = loadStep;

NRiter = 0;                                 % iteration count

while NRiter < beam.maxNumIter              % LOOP NR AT CURRENT LOAD STEP
    NRiter = NRiter + 1;

    [x, deltax] = solveLinSys(beam,mesh,alpha,x);

    if norm(deltax) <= beam.NRtol
        break
    end
end                                         % END LOOP NR AT CURRENT LOAD STEP

% storing NR-converged solutions at the current load step
mesh.num_iters(i)   = NRiter;
beam.q(:,i+1)       = x(1:mesh.ndofs);
beam.e(:,i+1)       = x(mesh.ndofs+1:mesh.ndofs+mesh.ndofse);
beam.s(:,i+1)       = x(mesh.ndofs+mesh.ndofse+1:mesh.ndofs+mesh.ndofse+mesh.ndofss);
beam.chi(:,i+1)     = x(mesh.ndofs+mesh.ndofse+mesh.ndofss+1:end);

% stopping computation in case of non-converging NR solution at the current load step
if NRiter == beam.maxNumIter
    warning('Newton-Raphson scheme does not converge at load step %.1d.\nComputation stopped.\n', i);
    beam.NRflag = 1;
end