function K = aseembleLinTangentKrollingM(beam, mesh, x, alpha, computeAnumerically)

% allocation
K = zeros(mesh.ndofs);

% assemble the tangent stiffness matrix of the rolling moment term
if computeAnumerically
    epsilon = 1e-8;
    for i = 1:mesh.ndofs
        x1 = x; x2 = x;
        x1(i) = x1(i) - epsilon;
        x2(i) = x2(i) + epsilon;

        f1 = assembleFextRollingM(beam, mesh, alpha, x1);
        f2 = assembleFextRollingM(beam, mesh, alpha, x2);

        K(:,i) = (f2 - f1) ./ (2*epsilon);
    end    
else
    % rollup moment at the current load step
    m = beam.MmaxEnd .* alpha;
    
    eleId = mesh.nelms;
    dofsInd = mesh.LM(:,eleId);

    % based on [Betsch, 2002, IJNME, Eq. 7]
    % jk-block: Kjk = mi (delta dj * Delta dk)
    for i = 1:3
        if i == 1
            rowIds = dofsInd(19:21);          % delta d2
            colIds = dofsInd(22:24);          % Delta d3
        elseif i == 2
            rowIds = dofsInd(22:24);          % delta d3
            colIds = dofsInd(16:18);          % Delta d1
        else
            rowIds = dofsInd(16:18);          % delta d1
            colIds = dofsInd(19:21);          % Delta d2
        end
        K(rowIds,colIds) = m(i) .* eye(3);
    end
    

    % based on [Gebhardt, 2020, Journal of Nonlinear Science]
%     ids1 = dofsInd(16:18);
%     ids2 = dofsInd(end-2:end);
%     K(ids1,ids2) = 0.5*m .* eye(3);
%     K(ids2,ids1) = -0.5*m .* eye(3);
end

K = sparse(K);