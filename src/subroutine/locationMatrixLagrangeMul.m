function mesh = locationMatrixLagrangeMul(mesh)

numLagrangeMulPerEle = mesh.numLagrangeMulPerNode * 2;
mesh.LMlamu = zeros(numLagrangeMulPerEle, mesh.nelms);

mesh.IEN = elementNodeMatrix(mesh);
mesh.IDlamu = destinationArrayLagrangeMul(mesh);

for e = 1:mesh.nelms
    n1 = mesh.IEN(1,e);
    n2 = mesh.IEN(2,e);
    mesh.LMlamu(:,e) = [mesh.IDlamu(:,n1); mesh.IDlamu(:,n2)];
end