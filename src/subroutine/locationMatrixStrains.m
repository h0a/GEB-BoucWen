function mesh = locationMatrixStrains(mesh)

numDofsPerEle = 6;      % 6 strain components
mesh.LMstrains = zeros(numDofsPerEle, mesh.nelms);

% shifting for the mixed formulation
n = mesh.ndofs;

for e = 1:mesh.nelms    
    mesh.LMstrains(:,e) = (1:numDofsPerEle) + (e-1)*numDofsPerEle + n;
end