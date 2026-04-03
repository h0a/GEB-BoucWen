function fext = assembleFextBody(beam, mesh)

fext = zeros(mesh.ndofs,1);

for e = 1:mesh.nelms
    dofsInd = mesh.LM(:,e);
    fext(dofsInd) = fext(dofsInd) + elementBodyForceVector(beam, mesh, e);
end