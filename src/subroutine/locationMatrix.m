function mesh = locationMatrix(mesh)

numDofsPerEle = mesh.numDofsPerNode * 2;
mesh.LM = zeros(numDofsPerEle, mesh.nelms);

mesh.IEN = elementNodeMatrix(mesh);
mesh.ID = destinationArray(mesh);

for e = 1:mesh.nelms
    n1 = mesh.IEN(1,e);
    n2 = mesh.IEN(2,e);
    mesh.LM(:,e) = [mesh.ID(:,n1); mesh.ID(:,n2)];
end