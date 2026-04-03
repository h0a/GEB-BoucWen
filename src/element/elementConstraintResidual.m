function r = elementConstraintResidual(mesh, eleId, x)
% compute the constraint residual
% input:  mesh:      struct incl. relevant information of the discretization
%         eleId:     index of the current element
%         x:         vector of the solution of the current configuration
% output: r:         constraint residual vector 12x1 (2x(6 constraint per node) per element)


% (nodal) q of the current element
qe = x(mesh.LM(:,eleId));

% q at left and right node
qleft = qe(1:mesh.numDofsPerNode);
qright = qe(mesh.numDofsPerNode+1:end);

% constraint vector at each node
hleft = computeHVector(qleft);
hright = computeHVector(qright);

% constraint residual - factor 0.5 due to repeating at each node during the assembly routine - except the first and last node
if mesh.nelms == 1
    r = [hleft; hright];
else
    if eleId == 1
        r = [hleft; 0.5 .* hright];
    elseif eleId == mesh.nelms
        r = [0.5 .* hleft; hright];
    else
        r = 0.5 .* [hleft; hright];
    end
end