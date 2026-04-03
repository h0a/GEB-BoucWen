function q = assembleGlobalQ(mesh, nodalphi0, nodald1, nodald2, nodald3)

% assemble global vector q of all nodes

q = zeros(mesh.ndofs,1);
for n = 1:mesh.numNodes
    id = (mesh.numDofsPerNode*(n-1)+1):(mesh.numDofsPerNode*n);
    q(id) = [nodalphi0(:,n); nodald1(:,n); nodald2(:,n); nodald3(:,n)];
end