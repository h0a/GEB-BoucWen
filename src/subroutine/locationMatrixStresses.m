function mesh = locationMatrixStresses(mesh)

numDofsPerEle = 6;      % 6 stress components
mesh.LMstresses = zeros(numDofsPerEle, mesh.nelms);

% shifting for the mixed formulation
n = mesh.ndofs + mesh.ndofse;

for e = 1:mesh.nelms    
    mesh.LMstresses(:,e) = (1:numDofsPerEle) + (e-1)*numDofsPerEle + n;
end