function zrelax = relaxBWz(beam, mesh, zPre, zCur)
% relax for each z component corresponding to each strain dof

% allocation
zrelax = zeros(numel(zCur),1);

% elementwise relaxation
for e = 1:mesh.nelms                          % LOOP OVER ELEMENTS
    % corresponding z of the current element
    zPe = zPre(mesh.LMstrains(:,e)-mesh.ndofs);
    zCe = zCur(mesh.LMstrains(:,e)-mesh.ndofs);

    % relaxation for z-components
    zrelax(mesh.LMstrains(:,e)-mesh.ndofs) = (1 - beam.BWomega) .* zPe + beam.BWomega .* zCe;
end                                           % END LOOP OVER ELEMENTS