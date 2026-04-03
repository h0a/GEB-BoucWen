function ID = destinationArrayLagrangeMul(mesh)

% ID(dc,n) = dg; 
% where: dc is the local dof number at the n-th node
%        dg is the global dof number

% standard numbering
ID = zeros(mesh.numLagrangeMulPerNode,mesh.numNodes);

for n = 1:mesh.numNodes
    ID(:,n) = (1:mesh.numLagrangeMulPerNode) + (n-1)*mesh.numLagrangeMulPerNode;
end

% shift to the last variable field of the mixed formulation
n = mesh.ndofs + mesh.ndofse + mesh.ndofss;
ID = ID + n;