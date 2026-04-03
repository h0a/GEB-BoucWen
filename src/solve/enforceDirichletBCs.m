function [A,rhs] = enforceDirichletBCs(A,rhs,mesh)

% collect all fixed dofs at both 2 boundaries
removedDofIds = [];
for i = 1:2
    removedDofIds = [removedDofIds;mesh.fixedDofsInds{i};mesh.removedLamuInds{i}];
end
removedDofIds = sort(removedDofIds);

% eliminate rows and columns of fixed dofs and irrelevant Lagrange multipliers
rhs(removedDofIds) = [];
A(removedDofIds,:) = [];
A(:,removedDofIds) = [];