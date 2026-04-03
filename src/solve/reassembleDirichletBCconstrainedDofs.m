function xFull = reassembleDirichletBCconstrainedDofs(xConstrained,mesh)

if length(xConstrained) == mesh.numActiveDofs    
    xFull = zeros(mesh.ndofsTotal,1);

    removedDofIds = [];
    for i = 1:2
        removedDofIds = [removedDofIds; mesh.fixedDofsInds{i};mesh.removedLamuInds{i}];
    end
   
    ids = sort(setdiff(1:mesh.ndofsTotal,removedDofIds));

    if length(ids) == length(xConstrained)
        xFull(ids) = xConstrained;
    else
        fprintf('\nUnmatching length of the solution and the number of active dofs. Return the constrained solution.\n');
        xFull = xConstrained;
    end
else
    fprintf('\nUnmatching length of the solution and the number of active dofs. Return the constrained solution.\n');
    xFull = xConstrained;
end
