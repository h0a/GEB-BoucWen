function K = aseembleLinTangentKfollowingF(beam, mesh, x, alpha, computeAnumerically)

% allocation
K = zeros(mesh.ndofs);

% assemble the tangent stiffness matrix of the rolling moment term
if computeAnumerically
    epsilon = 1e-8;
    for i = 1:mesh.ndofs
        x1 = x; x2 = x;
        x1(i) = x1(i) - epsilon;
        x2(i) = x2(i) + epsilon;

        f1 = assembleFextFollowingF(beam, mesh, alpha, x1);
        f2 = assembleFextFollowingF(beam, mesh, alpha, x2);

        K(:,i) = (f2 - f1) ./ (2*epsilon);
    end    
else
    % vector of the constant external moment (at the current load step)
    f = beam.Fmax * alpha;
    
    eleId = mesh.nelms;
    dofsInd = mesh.LM(:,eleId);
    ids1 = dofsInd(13:15);
    ids2 = dofsInd(16:18);

    K(ids1,ids2) = f .* eye(3);
end

K = sparse(K);