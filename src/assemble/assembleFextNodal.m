function fext = assembleFextNodal(beam, mesh)

fext = zeros(mesh.ndofs,1);

for n = 1:mesh.numNodes
    dofsInd = mesh.ID(:,n);
    fext(dofsInd(1:3)) = beam.nodalFext(:,n);
end