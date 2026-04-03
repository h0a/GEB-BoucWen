function [nodalphi0, nodald1, nodald2, nodald3] = deassembleGlobalQ(q, mesh)

% deassemble global vector q of all nodes to position and nodal directors 
% of the format (3 coordinates x number of nodes)

nodalphi0 = zeros(3,mesh.numNodes);
nodald1 = zeros(3,mesh.numNodes);
nodald2 = zeros(3,mesh.numNodes);
nodald3 = zeros(3,mesh.numNodes);

for n = 1:mesh.numNodes
    id = (mesh.numDofsPerNode*(n-1)+1):(mesh.numDofsPerNode*n);
    qn = q(id);
    nodalphi0(:,n) = qn(1:3);
    nodald1(:,n) = qn(4:6);
    nodald2(:,n) = qn(7:9);
    nodald3(:,n) = qn(10:12);
end