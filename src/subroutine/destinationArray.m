function ID = destinationArray(mesh)

% ID(dc,n) = dg; 
% where: dc is the local dof number at the n-th node
%        dg is the global dof number

% standard numbering
ID = zeros(mesh.numDofsPerNode,mesh.numNodes);

for n = 1:mesh.numNodes
    ID(:,n) = (1:mesh.numDofsPerNode) + (n-1)*mesh.numDofsPerNode;
end